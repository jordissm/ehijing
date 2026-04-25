#pragma once

#include <cstdint>
#include <string>

/**
 * @brief Configuration parameters for an eHIJING/ELECTRA run.
 *
 * This structure aggregates all runtime parameters required to generate
 * events and manage output in a reproducible and chunked workflow.
 *
 * It is typically populated via command-line parsing (see parse_args()) and
 * then passed to the main execution logic of the program.
 *
 * The configuration supports:
 * - Event generation control (number of events, seeds)
 * - Nuclear system definition (Z, A)
 * - Output organization (run directories, chunking)
 * - Physics configuration (TMD parameters, external config files)
 */
struct RunConfig {

    /// Total number of events to generate in this invocation.
    int n_events;

    /// Atomic number \f$ Z \f$ of the target nucleus.
    int atomic_number;

    /// Mass number \f$ A \f$ of the target nucleus.
    int mass_number;

    /// Constant factor controlling the transverse-momentum-dependent (TMD) scale.
    double tmd_k_constant;

    /// Path to precomputed tables.
    std::string tabulation_path;

    /// Root directory where event outputs will be written.
    ///
    /// Typically contains per-event or per-chunk subdirectories, along with
    /// metadata files and completion markers.
    std::string run_path;

    /// Path to the hard-process configuration file (e.g., PYTHIA settings).
    std::string hard_process_config_path;

    /// Path to the hadronization configuration file (e.g., fragmentation settings).
    std::string hadronization_config_path;

    /// Mode for medium modification of the parton shower.
    int medium_modification_mode;

    /// ID of the first event in this run.
    ///
    /// Used for deterministic seeding and for splitting large productions into
    /// independent chunks.
    int64_t first_event_id;

    /// Number of events per chunk.
    ///
    /// This controls how events are grouped into directories (e.g.,
    /// `events_00000000-00000999/`) and is important for parallel workflows.
    int64_t chunk_size;

    /// Base random seed for event generation.
    ///
    /// Typically combined with the event ID to ensure reproducibility across runs.
    uint32_t seed;
};

/**
 * @brief Print program usage information and exit.
 *
 * Displays a help message describing the expected command-line arguments and
 * their meaning, then terminates the program.
 *
 * @param prog Name of the executable (usually argv[0]).
 */
void usage(const char* prog);

/**
 * @brief Parse command-line arguments into a RunConfig structure.
 *
 * Interprets the provided command-line arguments and constructs a RunConfig
 * object containing all required runtime parameters.
 *
 * This function is responsible for:
 * - Validating required arguments
 * - Applying default values where appropriate
 * - Converting string inputs into typed configuration fields
 *
 * @param argc Argument count from main().
 * @param argv Argument vector from main().
 *
 * @return A fully populated RunConfig structure.
 *
 * @throws std::runtime_error if required arguments are missing or invalid.
 *
 * @note The exact set of supported command-line options should match the
 *       usage() output.
 */
RunConfig parse_args(int argc, char* argv[]);