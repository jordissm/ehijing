#include "cli.hpp"

#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>

namespace {

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
        << "  " << prog
        << " --nevents N --first-event-id I --Z Z --A A --mode M --K K "
        << "--table-dir PATH --run-dir PATH --config-file PATH --seed S\n\n"
        << "Example:\n"
        << "  " << prog
        << " --nevents 1000 --first-event-id 0 --Z 1 --A 2 --mode 0 --K 4.0 "
        << "--table-dir output/runs/ehijing/tables/K4p0 "
        << "--run-dir output/runs/ehijing/events "
        << "--config-file /opt/electra/ehijing_bin/config/experiments/hermes.setting "
        << "--seed 12345\n";
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

    // Required arguments parser
    cfg.nEvents    = std::stoi(require_arg(args, "--nevents", argv[0]));
    cfg.Z          = std::stoi(require_arg(args, "--Z", argv[0]));
    cfg.A          = std::stoi(require_arg(args, "--A", argv[0]));
    cfg.mode       = std::stoi(require_arg(args, "--mode", argv[0]));
    cfg.K          = std::stod(require_arg(args, "--K", argv[0]));
    cfg.tableDir   = require_arg(args, "--table-dir", argv[0]);
    cfg.runDir     = require_arg(args, "--run-dir", argv[0]);
    cfg.configFile = require_arg(args, "--config-file", argv[0]);

    // Parse optional first global event ID for this chunk
    cfg.firstEventId = 0;
    if (auto it = args.find("--first-event-id"); it != args.end()) {
        cfg.firstEventId = std::stoll(it->second);
        if (cfg.firstEventId < 0) {
            std::cerr << "ERROR: --first-event-id must be >= 0\n";
            std::exit(2);
        }
    }

    // Parse optional chunk size argument
    // By default, run all events in one chunk
    cfg.chunkSize = cfg.nEvents;
    if (auto it = args.find("--chunk-size"); it != args.end()) {
        cfg.chunkSize = std::stoll(it->second);
        if (cfg.chunkSize <= 0) {
            std::cerr << "ERROR: --chunk-size must be > 0\n";
            std::exit(2);
        }
    }

    // Parse optional random seed argument for reproducibility
    if (auto it = args.find("--seed"); it != args.end()) {
        // Pythia seeds must be in [1, 900000000] typically; keep it in range
        const uint64_t s64 = std::stoull(it->second);
        cfg.seed = static_cast<uint32_t>(1 + (s64 % 900000000ULL));
    } else {
        // Non-deterministic fallback
        cfg.seed = static_cast<uint32_t>(std::random_device{}());
    }

    if (cfg.nEvents <= 0) {
        std::cerr << "ERROR: --nevents must be > 0\n";
        std::exit(2);
    }

    if (cfg.A <= 0) {
        std::cerr << "ERROR: --A must be > 0\n";
        std::exit(2);
    }

    if (cfg.Z < 0 || cfg.Z > cfg.A) {
        std::cerr << "ERROR: require 0 <= Z <= A\n";
        std::exit(2);
    }

    return cfg;
}