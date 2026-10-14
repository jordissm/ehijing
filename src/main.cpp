#include "Pythia8/Pythia.h"
#include "cli.hpp"
#include "filesystem_utils.hpp"
#include "dis_kinematics.hpp"
#include "ehijing_constants.hpp"
#include "event_output.hpp"
#include "hadronizer.hpp"
#include "modified_ff.hpp"
#include "spacetime_history.hpp"

using namespace Pythia8;
namespace constants = ehijing::constants;
#include <algorithm>
#include <chrono>
#include <cstdlib> // for rand() and srand()
#include <ctime>   // for time()
#include <fstream>
#include <random>
#include <memory>
#include <iomanip>
#include <sstream>
#include <unistd.h>
#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_map>


std::random_device rrd;  // Non-deterministic generator
std::mt19937 gen(rrd()); // Mersenne Twister engine

// Helper function to set Pythia settings
template <class T> void add_arg(Pythia& pythia, std::string name, T value) {
    std::stringstream ss;
    ss << name << " = " << value;
    std::cout << ss.str() << std::endl;
    pythia.readString(ss.str());
}


int main(int argc, char* argv[]) {

    const RunConfig cfg = parse_args(argc, argv);

    const int n_events                              = cfg.n_events;
    const int atomic_number                         = cfg.atomic_number;
    const int mass_number                           = cfg.mass_number;
    const double tmd_k_constant                     = cfg.tmd_k_constant;
    const int medium_modification_mode              = cfg.medium_modification_mode;
    const std::string& tabulation_path              = cfg.tabulation_path;
    const std::string& run_path                     = cfg.run_path;
    const std::string& hard_process_config_path     = cfg.hard_process_config_path;
    const std::string& hadronization_config_path    = cfg.hadronization_config_path;
    const std::string& dis_cuts_config_path         = cfg.dis_cuts_config_path;
    const int64_t first_event_id                    = cfg.first_event_id;
    const int64_t chunk_size                        = cfg.chunk_size;
    const uint32_t seed                             = cfg.seed;

    if (cfg.nuclear_target) {
        try {
            EHIJING::NuclearTarget validate(mass_number, atomic_number, cfg.nucleon_width_fm,
                cfg.nuclear_radius_fm, cfg.nuclear_surface_fm);
            if (EHIJING::nuclear_target_abi_version() != 1) return 2;
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << '\n'; return 2;
        }
    }
    // Ensure output directories exist
    try { // Attempt to create the directories if they do not exist
        std::filesystem::create_directories(std::filesystem::path(run_path));
        std::filesystem::create_directories(std::filesystem::path(tabulation_path));
    } catch (const std::filesystem::filesystem_error& e) { // Handle any errors that occur during directory creation
        std::cerr << "ERROR: cannot create output/table directories:\n  " << e.what() << "\n";
        return 3;
    }

    // Build the nucleus ID used by Pythia
    const int nucleus_PDG_id =
        constants::pythia::nucleus_pdg_base +
        atomic_number * constants::pythia::nucleus_pdg_z_factor +
        mass_number * constants::pythia::nucleus_pdg_a_factor;

    /*
    The nuclear modication to be used for beam A
        nuclear_modification = 0: Only Isospin effect.
        nuclear_modification = 1: EPS09, LO
        nuclear_modification = 2: EPS09, NLO
        nuclear_modification = 3: EPPS16, NLO
    */
    // We will use only isospin for deuteron, and EPPS16 NLO for heavier nucleus
    int nuclear_modification = (atomic_number > 2) ? 3 : 0;

    // Initialize the Pythia instance for hadronization
    Hadronizer hadronizer(hadronization_config_path);

    const DISCuts dis_cuts = load_dis_cuts(dis_cuts_config_path);

    // Create generator object for the eHIJING-Pythia high-Q parton shower in a medium
    Pythia pythia;

    // Read settings from the config file and apply them to the Pythia instance
    pythia.readFile(hard_process_config_path);
    
    // Make Pythia deterministic
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));

    // Override the default proton mass with the PDG value to ensure consistency with the nuclear geometry and kinematics
    pythia.particleData.m0(2212, 0.938);

    gen.seed(seed);
    std::srand(seed);

    // Set Pythia and eHIJING settings from command line arguments
    // Note: Command line arguments will override settings in the config file if there are conflicts
    add_arg<int>(pythia, "PDF:nPDFSetA", nuclear_modification);
    add_arg<int>(pythia, "PDF:nPDFBeamA", nucleus_PDG_id);
    add_arg<int>(pythia, "eHIJING:MediumModificationSplittingFunctionMode", medium_modification_mode);
    add_arg<int>(pythia, "eHIJING:MassNumber", mass_number);
    add_arg<int>(pythia, "eHIJING:AtomicNumber", atomic_number);
    add_arg<double>(pythia, "eHIJING:TransverseMomentumDistributionConstant", tmd_k_constant);
    add_arg<std::string>(pythia, "eHIJING:TabulationPath", tabulation_path);

    add_arg<std::string>(pythia, "eHIJING:ExplicitNucleus", cfg.nuclear_target ? "on" : "off");
    add_arg<double>(pythia, "eHIJING:NucleonWidth", cfg.nucleon_width_fm);
    add_arg<double>(pythia, "eHIJING:NuclearRadius", cfg.nuclear_radius_fm);
    add_arg<double>(pythia, "eHIJING:NuclearSurface", cfg.nuclear_surface_fm);
    // Initialize Pythia
    if (!pythia.init()) {
        std::cerr << "ERROR: PYTHIA initialization failed\n"; return 2;
    }

    // Initialize the modified FF module for medium corrections to the parton shower
    double xg_n = pythia.settings.parm("eHIJING:xG-n");
    double xg_lambda = pythia.settings.parm("eHIJING:xG-lambda");
    Modified_FF modified_ff(medium_modification_mode, 
                            atomic_number, 
                            mass_number, 
                            tmd_k_constant, 
                            xg_n,
                            xg_lambda, 
                            tabulation_path);

    // Define counter for events written to output
    int n_written = 0;

    // Define counter for events that failed to generate
    int n_generation_failed = 0;
    int n_recoil_rejected = 0, n_balance_rejected = 0;
    int attempts_since_written = 0;

    // Define counter for events that failed the DIS trigger
    int n_trigger_failed = 0;

    // Define counter for events that failed hadronization
    const int hadronization_retries_max = 15; // Maximum number of hadronization failures before giving up
    int n_hadronization_failed = 0;

    // Define counter for total events
    int n_total_attempts = 0;
    
    // Begin event loop
    while (n_written < n_events) {

        // Flag to indicate whether the event was successfully generated, passed the trigger, and hadronized
        bool success = false;

        for (int attempt = 0; attempt < hadronization_retries_max; ++attempt) { // Try to generate a valid event

            // Add to total event attempt count
            n_total_attempts++;
            if (++attempts_since_written > 10000 && cfg.nuclear_target) {
                std::cerr << "ERROR: 10000 consecutive attempts produced no accepted event; "
                    << "generation_failures=" << n_generation_failed << " recoil_rejections=" << n_recoil_rejected
                    << " balance_rejections=" << n_balance_rejected << " trigger_failures=" << n_trigger_failed
                    << " hadronization_failures=" << n_hadronization_failed << '\n';
                return 4;
            }

            // Skip event if it fails to generate
            if (!pythia.next()) {
                // Add to failed event generation count
                n_generation_failed++;
                continue;
            };

            if (cfg.nuclear_target && !pythia.event.nuclearTarget()) {
                std::cerr << "ERROR: explicit nuclear target missing from generated event; "
                          << "check rebuilt PYTHIA library/settings and supported DIS configuration\n";
                return 2;
            }
            // Skip event if its kinematics do not satisfy the DIS trigger conditions
            const DISKinematics kinematics = compute_dis_kinematics(pythia);
            if (!is_valid_dis_event(kinematics, dis_cuts)) {
                // Add to failed DIS trigger count
                n_trigger_failed++;
                continue;
            }

            // Find the hard vertex
            double Rx = pythia.event.Rx();
            double Ry = pythia.event.Ry();
            double Rz = pythia.event.Rz();
            std::unique_ptr<SpacetimeHistory> history;
            if (cfg.history) {
                history = std::make_unique<SpacetimeHistory>(first_event_id + n_written, attempt + 1);
                history->set_context(Vec4(Rx, Ry, Rz, 0.0), kinematics.pProton,
                                     pythia.settings.flag("PartonVertex:setVertex"));
                history->capture(pythia.process, "hard_process", "process", "source_parton_vertices");
                history->capture(pythia.event, "high_scale_shower", "shower", "source_parton_vertices");
            }
            if (history && pythia.event.nuclearTarget()) {
                const auto target = pythia.event.nuclearTarget();
                std::vector<Particle> nucleons;
                for (size_t id = 0; id < target->nucleons().size(); ++id) {
                    const auto& n = target->nucleons()[id];
                    const double mass = n.proton ? constants::mass::proton_gev : constants::mass::neutron_gev;
                    Vec4 p(0.,0.,0.,mass); p.bst(kinematics.pProton);
                    Particle observed(n.proton ? 2212 : 2112,2,0,0,0,0,0,0,p,mass);
                    Vec4 v(n.position[0],n.position[1],n.position[2],0.);
                    v *= constants::units::hbarc_gev_fm*constants::units::fm_to_mm;
                    v.bst(kinematics.pProton); observed.vProd(v); observed.nuclearTargetId(id);
                    nucleons.push_back(observed);
                }
                history->capture_particles(nucleons,"nuclear_target","nucleus","absolute_mm_nuclear_origin",{}, {},
                    "Frozen Gaussian nucleons; width_fm=" + std::to_string(target->width()*constants::units::hbarc_gev_fm)
                    + "; primary_target_id=" + std::to_string(target->primary()));
            }
            // Modify the final shower with low-Q^2 medium corrections
            if (!modified_ff.sample_ff_partons(pythia, kinematics, Rx, Ry, Rz, history.get())) {
                ++n_recoil_rejected; continue;
            }
            if (history) history->capture(pythia.event, "medium_modified", "shower",
                                          "relative_mm_hard_interaction");
            
            // Put the parton-level event into the separate hadronizer
            auto hadronized_event_opt = hadronizer.hadronize(
                pythia, atomic_number, mass_number, Rx, Ry, Rz,
                first_event_id + n_written, attempt + 1, history.get());
            if (!hadronized_event_opt) {
                std::cerr << "WARNING: hadronization failed for event " << (first_event_id + n_written) << " on attempt " << (attempt + 1) << "/" << hadronization_retries_max << "\n";
                n_hadronization_failed++;
                continue;
            }
            std::cout << "INFO: hadronization succeeded for event " << (first_event_id + n_written) << " on attempt " << (attempt + 1) << "/" << hadronization_retries_max << "\n";
            const auto& hadronization_result = *hadronized_event_opt;
            
            const auto target = pythia.event.nuclearTarget();
            const auto balance = calculate_spectator_balance(atomic_number,mass_number,
                hadronization_result.struck_nucleon,hadronization_result.particles);
            if (target && !target_balance_matches(*target,balance)) {
                std::cerr << "WARNING: rejecting event: target ledger disagrees with final baryon number/charge\n";
                ++n_balance_rejected; continue;
            }
            // Set event ID
            const int64_t event_id = first_event_id + n_written;
            
            // All event outputs share the same event ID and shard directory.
            const EventPaths paths = make_event_paths(std::filesystem::path(run_path), event_id, chunk_size);

            // Open OSCAR output file
            std::ofstream event_out(paths.eventPath);
            if (!event_out) {
                std::cerr << "ERROR: cannot open output file: " << paths.eventPath << std::endl;
                return 1;
            }
            
            // Open metadata output file
            std::ofstream meta_out(paths.metaPath);
            if (!meta_out) {
                std::cerr << "ERROR: cannot open metadata file: " << paths.metaPath << std::endl;
                return 1;
            }

            // Only the successful hadronization attempt is paired with this OSCAR event.
            if (!hadronization_result.fragmentation_trace.empty()) {
                std::ofstream fragmentation_out(paths.fragmentationPath);
                fragmentation_out << hadronization_result.fragmentation_trace;
                fragmentation_out.close();
                if (!fragmentation_out) {
                    std::cerr << "ERROR: cannot write fragmentation file: "
                              << paths.fragmentationPath << std::endl;
                    return 1;
                }
            }

            // Write the event data and metadata to the respective files
            write_event_headers(event_out);
            write_event_output(event_id,
                               atomic_number,
                               mass_number,
                               hadronization_result.struck_nucleon,
                               kinematics,
                               hadronization_result.particles,
                               event_out,
                               meta_out,
                               history.get(), target.get());

            if (history) {
                try {
                    std::ofstream history_out(paths.historyPath);
                    if (!history_out) throw std::runtime_error("Cannot open history output");
                    history->write(history_out);
                    history_out.close();
                    if (!history_out) throw std::runtime_error("Cannot close history output");
                } catch (const std::exception& error) {
                    std::cerr << "ERROR: " << error.what() << ": " << paths.historyPath << '\n';
                    return 1;
                }
            }

            // Mark event as successfully generated, passed trigger, and hadronized
            success = true;
            n_written++;
            attempts_since_written = 0;
            break; // Exit the retry loop if the event was successful
        }

    }

    std::cout << "INFO: written=" << n_written << " attempts=" << n_total_attempts
        << " recoil_rejections=" << n_recoil_rejected << " balance_rejections=" << n_balance_rejected
        << " generation_failures=" << n_generation_failed << " trigger_failures=" << n_trigger_failed
        << " hadronization_failures=" << n_hadronization_failed << '\n';
    // Done
    return 0;
}
