#include "Pythia8/Pythia.h"
#include "cli.hpp"
#include "filesystem_utils.hpp"
#include "dis_kinematics.hpp"
#include "event_output.hpp"
#include "hadronizer.hpp"
#include "modified_ff.hpp"

using namespace Pythia8;
#include <algorithm>
#include <chrono>
#include <cstdlib> // for rand() and srand()
#include <ctime>   // for time()
#include <fstream>
#include <random>
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
    const int64_t first_event_id                    = cfg.first_event_id;
    const int64_t chunk_size                        = cfg.chunk_size;
    const uint32_t seed                             = cfg.seed;

    // Ensure output directories exist
    try { // Attempt to create the directories if they do not exist
        std::filesystem::create_directories(std::filesystem::path(run_path));
        std::filesystem::create_directories(std::filesystem::path(tabulation_path));
    } catch (const std::filesystem::filesystem_error& e) { // Handle any errors that occur during directory creation
        std::cerr << "ERROR: cannot create output/table directories:\n  " << e.what() << "\n";
        return 3;
    }

    // Build the nucleus ID used by Pythia
    const int nucleus_PDG_id = 100000000 + atomic_number * 10000 + mass_number * 10;

    /*
    The nuclear modication to be used for beam A
        nuclear_modification = 0: Only Isospin effect.
        nuclear_modification = 1: EPS09, LO
        nuclear_modification = 2: EPS09, NLO
        nuclear_modification = 2: EPPS16, NLO
    */
    // We will use only isospin for deuteron, and EPPS16 NLO for heavier nucleus
    int nuclear_modification = (atomic_number > 2) ? 3 : 0;

    // Initialize the Pythia instance for hadronization
    Hadronizer hadronizer;

    // Create generator object for the eHIJING-Pythia high-Q parton shower in a medium
    Pythia pythia;

    // Read settings from the config file and apply them to the Pythia instance
    pythia.readFile(hard_process_config_path);
    
    // Make Pythia deterministic
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));

    gen.seed(seed);
    std::srand(seed);

    // Set Pythia and eHIJING settings from command line arguments
    // Note: Command line arguments will override settings in the config file if there are conflicts
    add_arg<int>(pythia, "PDF:nPDFSetA", nuclear_modification);
    add_arg<int>(pythia, "PDF:nPDFBeamA", nucleus_PDG_id);
    add_arg<int>(pythia, "eHIJING:MediumModificationSplittingFunctionMode", medium_modification_mode);
    add_arg<int>(pythia, "eHIJING:MassNumber", mass_number);
    add_arg<int>(pythia, "eHIJING:ChargeNumber", atomic_number);
    add_arg<double>(pythia, "eHIJING:Kfactor", tmd_k_constant);
    add_arg<std::string>(pythia, "eHIJING:TablePath", tabulation_path);

    // Initialize Pythia
    pythia.init();

    // Initialize the modified FF module for medium corrections to the parton shower
    double xg_n = pythia.settings.parm("eHIJING:xG-n");
    double xg_lambda = pythia.settings.parm("eHIJING:xG-lambda");
    Modified_FF modified_ff(mode, 
                            atomic_number, 
                            mass_number, 
                            k_factor, 
                            xg_n,
                            xg_lambda, 
                            table_path);

    // Define counter for events written to output
    int n_written = 0;

    // Define counter for events that failed to generate
    int n_generation_failed = 0;

    // Define counter for events that failed the DIS trigger
    int n_trigger_failed = 0;

    // Define counter for events that failed hadronization
    const int hadronization_retries_max = 10; // Maximum number of hadronization failures before giving up
    int n_hadronization_failed = 0;

    // Define counter for total events
    int n_total_attempts = 0;
    
    // Begin event loop
    while (n_written < n_events) {

        // Flag to indicate whether the event was successfully generated, passed the trigger, and hadronized
        bool success = false;

        for (int attempt = 0; attempt < hadronization_retries_max; ++attempt) { // Try up to 10 times to generate a valid event

            // Add to total event attempt count
            n_total_attempts++;

            // Skip event if it fails to generate
            if (!pythia.next()) {
                // Add to failed event generation count
                n_generation_failed++;
                continue;
            };

            // Skip event if its kinematics do not satisfy the DIS trigger conditions
            const DISKinematics kinematics = compute_dis_kinematics(pythia);
            if (!is_valid_dis_event(kinematics)) {
                // Add to failed DIS trigger count
                n_trigger_failed++;
                continue;
            }

            // Find the hard vertex
            double Rx, Ry, Rz;
            // Modify the final shower with low-Q^2 medium corrections
            modified_ff.sample_ff_partons(pythia, Rx, Ry, Rz);
            
            // Put the parton-level event into the separate hadronizer
            auto hadronized_event_opt = hadronizer.hadronize(pythia, atomic_number, mass_number, Rx, Ry, Rz);
            if (!hadronized_event_opt) {
                n_hadronization_failed++;
                continue;
            }
            const auto& hadronized_event = *hadronized_event_opt;
            
            // Set event ID
            const int64_t event_id = first_event_id + n_written;
            
            // Define output paths for the OSCAR event file and the metadata file based on the event ID and chunk size
            const EventPaths paths = make_event_paths(std::filesystem::path(outDir), event_id, chunk_size);

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

            // Write the event data and metadata to the respective files
            write_event_headers(event_out);
            write_event_output(event_id, atomic_number, mass_number, kinematics, hadronized_event, event_out, meta_out);

            // Mark event as successfully generated, passed trigger, and hadronized
            success = true;
            n_written++;
            break; // Exit the retry loop if the event was successful
        }

    }

    // Done
    return 0;
}
