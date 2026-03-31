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

// Define a distribution (range 0.0 to 1.0)
std::uniform_real_distribution<> Ran_gen(0.0, 1.0);


double fermi_distribution(double r, double R_WS, double a_WS) {
    double f = 1. / (1. + exp((r - R_WS) / a_WS));
    return (f);
}

// S. Chakraborty et al., Physical Review C 107, 064318 (2023)
double sample_r_from_woods_saxon(double R_WS, double a_WS) {
    double rmaxCut = R_WS + 10. * a_WS;
    double r = 0.;
    do
    {
        r = rmaxCut * pow(Ran_gen(gen), 1.0 / 3.0);
    } while (Ran_gen(gen) > fermi_distribution(r, R_WS, a_WS));
    return (r);
}


// Helper function to set Pythia settings
template <class T> void add_arg(Pythia& pythia, std::string name, T value) {
    std::stringstream ss;
    ss << name << " = " << value;
    std::cout << ss.str() << std::endl;
    pythia.readString(ss.str());
}


int main(int argc, char* argv[]) {

    const RunConfig cfg = parse_args(argc, argv);

    const int nEvents = cfg.nEvents;
    const int Z = cfg.Z;
    const int A = cfg.A;
    const int mode = cfg.mode;
    const double K = cfg.K;
    const std::string& tableDir = cfg.tableDir;
    const std::string& outDir = cfg.runDir;
    const std::string& configFile = cfg.configFile;
    const int64_t first_event_id = cfg.firstEventId;
    const int64_t chunk_size = cfg.chunkSize;
    const uint32_t seed = cfg.seed;

    // Ensure output directories exist
    try { // Attempt to create the directories if they do not exist
        std::filesystem::create_directories(std::filesystem::path(outDir));
        std::filesystem::create_directories(std::filesystem::path(tableDir));
    } catch (const std::filesystem::filesystem_error& e) { // Handle any errors that occur during directory creation
        std::cerr << "ERROR: cannot create output/table directories:\n  " << e.what() << "\n";
        return 3;
    }

    // Build the nucleus ID used by Pythia & the PDF (the isospin effect)
    const int iNuclei = 100000000 + Z * 10000 + A * 10;

    // JORDI: These lines are different from ehijing-default-Briet-frame.cpp
    // Shadowing effect:
    //   nPDFset=0: only isospin
    //   nPDFset = 1: EPS09 LO
    //   nPDFset = 2: EPS09 NLO
    //   nPDFset = 2: EPPS16 NLO
    // We will use only isospin for deuteron,
    // and EPPS16 NLO for heavier nucleus
    int nPDFset = 0; // (A>2)?3:0;

    // Initialize the Pythia instance for hadronization
    Hadronizer hadronizer;

    // Initialize the eHIJING-Pythia for high-Q parton shower in medium
    Pythia pythia;

    // Read settings from the config file and apply them to the Pythia instance
    pythia.readFile(configFile);
    
    // JORDI: These lines are different from ehijing-default-Briet-frame.cpp
    // Make Pythia deterministic
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));

    // ALSO seed your own RNGs deterministically (important!)
    gen.seed(seed);
    std::srand(seed);

    // Set Pythia and eHIJING settings from command line arguments
    // Note: Command line arguments will override settings in the config file if there are conflicts
    add_arg<int>(pythia, "PDF:nPDFSetA", nPDFset);
    add_arg<int>(pythia, "PDF:nPDFBeamA", iNuclei);
    add_arg<int>(pythia, "eHIJING:Mode", mode);
    add_arg<int>(pythia, "eHIJING:AtomicNumber", A);
    add_arg<int>(pythia, "eHIJING:ChargeNumber", Z);
    add_arg<double>(pythia, "eHIJING:Kfactor", K);
    add_arg<std::string>(pythia, "eHIJING:TablePath", tableDir);

    // Prepare Pythia for event generation
    pythia.init();

    // Initialize the modified FF module for medium corrections to the parton shower
    Modified_FF MFF(mode, Z, A, K, pythia.settings.parm("eHIJING:xG-n"),
                    pythia.settings.parm("eHIJING:xG-lambda"), tableDir);

    // Define counter for triggered events
    int nTriggered = 0;

    // Define counter for failed events
    int nFailed = 0;

    // Define counter for total events
    int nTotal = 0;
    
    // Begin event loop
    while (nTriggered < nEvents) {

        // Add to total event count
        nTotal++;

        // Skip event if it fails to generate
        if (!pythia.next()) {
            // Add to failed event count
            nFailed++;
            continue;
        };

        // Skip event if its kinematics do not satisfy the DIS trigger conditions
        const DISKinematics kinematics = compute_dis_kinematics(pythia);
        if (!is_valid_dis_event(kinematics)) {
            continue;
        }

        // Event passed the trigger, add to triggered event count
        nTriggered++;

        // Set event ID
        const int64_t event_id = first_event_id + (nTriggered - 1);
        
        // Define output paths for the OSCAR event file and the metadata file based on the event ID and chunk size
        const EventPaths paths = make_event_paths(std::filesystem::path(outDir), event_id, chunk_size);
        
        // Initialize the (x, y, z) size
        double Rx, Ry, Rz;
        
        // Modify the final shower with low-Q^2 medium corrections
        MFF.sample_FF_partons(pythia.event, Rx, Ry, Rz);
        
        // Put the parton-level event into the separate hadronizer
        auto hadronizerEvent = hadronizer.hadronize(pythia, Z, A, Rx, Ry, Rz);
        
        // Open OSCAR output file
        std::ofstream fout_event(paths.eventPath);
        if (!fout_event) {
            std::cerr << "ERROR: cannot open output file: " << paths.eventPath << std::endl;
            return 1;
        }
        
        // Open metadata output file
        std::ofstream fout_meta(paths.metaPath);
        if (!fout_meta) {
            std::cerr << "ERROR: cannot open metadata file: " << paths.metaPath << std::endl;
            return 1;
        }

        // Write the event data and metadata to the respective files
        write_event_headers(fout_event);
        write_event_output(event_id, Z, A, kinematics, hadronizerEvent, fout_event, fout_meta);

        // Close the output files
        fout_event.close();
        fout_meta.close();

    }

    // Done
    return 0;
}