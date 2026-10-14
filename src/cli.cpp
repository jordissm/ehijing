#include "cli.hpp"
#include "ehijing_constants.hpp"

#include <cstdlib>
#include <charconv>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>

namespace {

namespace constants = ehijing::constants;

std::string require_arg(
    const std::unordered_map<std::string, std::string>& args,
    const std::string& key,
    const char* prog)
{
    auto it = args.find(key);
    if (it == args.end() || it->second.empty()) {
        std::cerr << "ERROR: missing required argument: " << key << "\n";
        usage(prog);
        std::exit(2);
    }
    return it->second;
}

} // namespace

void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " "
        << "--number-of-events <NUMBER OF EVENTS> "
        << "--first-event-id <ID> "
        << "--Z <ATOMIC NUMBER> "
        << "--A <MASS NUMBER> "
        << "--K <TMD CONSTANT> "
        << "--medium-modification-mode <MEDIUM MODIFICATION MODE (0: HT | 1: GHT)> "
        << "--tabulation-path <PATH> "
        << "--run-path <PATH> "
        << "--hard-process-config <PATH> "
        << "--hadronization-config <PATH> "
        << "--dis-cuts-config <PATH> "
        << "--seed <SEED> [--history <on|off>] [--nuclear-target <on|off>] [--nucleon-width-fm <fm>] [--nuclear-radius-fm <fm>] [--nuclear-surface-fm <fm>]\n\n";
}

RunConfig parse_args(int argc, char* argv[]) {
    std::unordered_map<std::string, std::string> args;

    for (int i = 1; i < argc; ++i) {
        std::string key = argv[i];

        if (key.rfind("--", 0) != 0) {
            std::cerr << "ERROR: unexpected positional arg: " << key << "\n";
            usage(argv[0]);
            std::exit(2);
        }

        if (i + 1 >= argc) {
            std::cerr << "ERROR: flag " << key << " requires a value\n";
            usage(argv[0]);
            std::exit(2);
        }

        args[key] = argv[++i];
    }

    RunConfig cfg{};

    if (auto it = args.find("--history"); it != args.end()) {
        if (it->second != "on" && it->second != "off") {
            std::cerr << "ERROR: --history must be on or off\n";
            std::exit(2);
        }
        cfg.history = it->second == "on";
    }

    if (auto it = args.find("--nuclear-target"); it != args.end()) {
        if (it->second != "on" && it->second != "off") {
            std::cerr << "ERROR: --nuclear-target must be on or off\n"; std::exit(2);
        }
        cfg.nuclear_target = it->second == "on";
    }
    if (args.count("--nucleon-width-fm")) cfg.nucleon_width_fm = std::stod(args.at("--nucleon-width-fm"));
    if (args.count("--nuclear-radius-fm")) cfg.nuclear_radius_fm = std::stod(args.at("--nuclear-radius-fm"));
    if (args.count("--nuclear-surface-fm")) cfg.nuclear_surface_fm = std::stod(args.at("--nuclear-surface-fm"));
    // Required arguments parser
    cfg.n_events        = std::stoi(require_arg(args, "--number-of-events", argv[0]));
    cfg.atomic_number   = std::stoi(require_arg(args, "--Z", argv[0]));
    cfg.mass_number     = std::stoi(require_arg(args, "--A", argv[0]));
    cfg.tmd_k_constant  = std::stod(require_arg(args, "--K", argv[0]));
    cfg.tabulation_path = require_arg(args, "--tabulation-path", argv[0]);
    cfg.run_path        = require_arg(args, "--run-path", argv[0]);
    cfg.hard_process_config_path    = require_arg(args, "--hard-process-config", argv[0]);
    cfg.hadronization_config_path   = require_arg(args, "--hadronization-config", argv[0]);
    cfg.dis_cuts_config_path        = require_arg(args, "--dis-cuts-config", argv[0]);
    cfg.medium_modification_mode    = std::stoi(require_arg(args, "--medium-modification-mode", argv[0]));

    // Parse optional first global event ID for this chunk
    cfg.first_event_id = 0;
    if (auto it = args.find("--first-event-id"); it != args.end()) {
        cfg.first_event_id = std::stoll(it->second);
        if (cfg.first_event_id < 0) {
            std::cerr << "ERROR: --first-event-id must be >= 0\n";
            std::exit(2);
        }
    }

    // Parse optional chunk size argument
    // By default, run all events in one chunk
    cfg.chunk_size = cfg.n_events;
    if (auto it = args.find("--chunk-size"); it != args.end()) {
        cfg.chunk_size = std::stoll(it->second);
        if (cfg.chunk_size <= 0) {
            std::cerr << "ERROR: --chunk-size must be > 0\n";
            std::exit(2);
        }
    }

    // Parse optional random seed argument for reproducibility
    if (auto it = args.find("--seed"); it != args.end()) {
        const auto& value = it->second;
        const auto parsed = std::from_chars(value.data(), value.data() + value.size(), cfg.seed);
        if (value.empty() || parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size()) {
            std::cerr << "ERROR: --seed must be an unsigned 64-bit decimal integer\n";
            std::exit(2);
        }
    } else {
        // Non-deterministic fallback
        std::random_device entropy;
        cfg.seed = (static_cast<uint64_t>(entropy()) << 32) | entropy();
    }

    if (cfg.n_events <= 0) {
        std::cerr << "ERROR: --number-of-events must be > 0\n";
        std::exit(2);
    }

    if (cfg.mass_number <= 0) {
        std::cerr << "ERROR: --A must be > 0\n";
        std::exit(2);
    }

    if (cfg.atomic_number < 0 || cfg.atomic_number > cfg.mass_number) {
        std::cerr << "ERROR: require 0 <= Z <= A\n";
        std::exit(2);
    }

    return cfg;
}
