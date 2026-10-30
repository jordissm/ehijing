#include "cli.hpp"
#include "ehijing_constants.hpp"
#include "nucleus_database.hpp"

#include <cstdlib>
#include <charconv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
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
        << "--nucleus <NAME> [--nucleus-database <YAML>] "
        << "--K <TMD CONSTANT> "
        << "--medium-modification-mode <MEDIUM MODIFICATION MODE (0: HT | 1: GHT)> "
        << "--tabulation-path <PATH> "
        << "--run-path <PATH> "
        << "--hard-process-config <PATH> "
        << "--hadronization-config <PATH> "
        << "--dis-cuts-config <PATH> "
        << "--seed <SEED> [--history <on|off>] [--nuclear-target <on|off>] [--nucleon-width-fm <fm>] [--nuclear-radius-fm <fm>] [--nuclear-surface-fm <fm>]\n"
        << "  " << prog << " --nucleus <NAME> [--nucleus-database <YAML>] --describe-nucleus\n"
        << "Legacy target arguments: --Z <Z> --A <A> instead of --nucleus.\n\n";
}

RunConfig parse_args(int argc, char* argv[]) {
    std::unordered_map<std::string, std::string> args;

    for (int i = 1; i < argc; ++i) {
        std::string key = argv[i];

        if (key == "--help") { usage(argv[0]); std::exit(0); }
        if (key == "--describe-nucleus") { args[key] = "on"; continue; }

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
    cfg.describe_nucleus = args.count("--describe-nucleus");
    if (args.count("--nucleus")) {
        if (args.count("--Z") || args.count("--A"))
            throw std::runtime_error("Use --nucleus or legacy --Z/--A, not both");
        cfg.nucleus = require_arg(args, "--nucleus", argv[0]);
        cfg.nucleus_database = ehijing::resolve_nucleus_database(
            args.count("--nucleus-database") ? require_arg(args, "--nucleus-database", argv[0]) : "",
            argv[0]).string();
        const auto target = ehijing::load_nucleus(cfg.nucleus_database, cfg.nucleus);
        cfg.atomic_number = target.Z; cfg.mass_number = target.A;
        cfg.distribution = target.distribution;
        cfg.nucleon_width_fm = target.nucleon_width_fm;
        cfg.nuclear_radius_fm = target.nuclear_radius_fm;
        cfg.nuclear_surface_fm = target.nuclear_surface_fm;
        cfg.hulthen_alpha_fm_inv = target.hulthen_alpha_fm_inv;
        cfg.hulthen_beta_fm_inv = target.hulthen_beta_fm_inv;
    } else {
        if (cfg.describe_nucleus || args.count("--nucleus-database"))
            throw std::runtime_error("--describe-nucleus and --nucleus-database require --nucleus");
        if (!args.count("--Z") && !args.count("--A"))
            throw std::runtime_error("Select a target with --nucleus NAME (or legacy --Z and --A)");
        cfg.atomic_number = std::stoi(require_arg(args, "--Z", argv[0]));
        cfg.mass_number = std::stoi(require_arg(args, "--A", argv[0]));
    }

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
    if (!(cfg.nucleon_width_fm > 0.) || !std::isfinite(cfg.nucleon_width_fm))
        throw std::runtime_error("--nucleon-width-fm must be finite and positive");
    if (!cfg.nucleus.empty()) {
        if (cfg.distribution == "hulthen" &&
            (args.count("--nuclear-radius-fm") || args.count("--nuclear-surface-fm")))
            throw std::runtime_error("Woods-Saxon overrides do not apply to a Hulthen target");
        if (cfg.distribution == "woods_saxon" &&
            (!(cfg.nuclear_radius_fm > 0.) || !std::isfinite(cfg.nuclear_radius_fm) ||
             !(cfg.nuclear_surface_fm > 0.) || !std::isfinite(cfg.nuclear_surface_fm)))
            throw std::runtime_error("Woods-Saxon radius and surface must be finite and positive");
    }
    if (cfg.describe_nucleus) return cfg;
    // Required arguments parser
    cfg.n_events        = std::stoi(require_arg(args, "--number-of-events", argv[0]));
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

std::string describe_nucleus(const RunConfig& cfg) {
    auto quote = [](const std::string& value) {
        std::ostringstream out;
        out << '"';
        for (unsigned char c : value) {
            if (c == '"' || c == '\\') out << '\\' << c;
            else if (c < 0x20) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c);
            else out << c;
        }
        out << '"';
        return out.str();
    };
    std::ostringstream out;
    out << std::setprecision(17)
        << "{\"nucleus\":" << quote(cfg.nucleus)
        << ",\"nucleus_database\":" << quote(cfg.nucleus_database)
        << ",\"distribution\":" << quote(cfg.distribution)
        << ",\"A\":" << cfg.mass_number << ",\"Z\":" << cfg.atomic_number
        << ",\"nucleon_width_fm\":" << cfg.nucleon_width_fm;
    if (cfg.distribution == "hulthen")
        out << ",\"hulthen_alpha_fm_inv\":" << cfg.hulthen_alpha_fm_inv
            << ",\"hulthen_beta_fm_inv\":" << cfg.hulthen_beta_fm_inv;
    else
        out << ",\"nuclear_radius_fm\":" << cfg.nuclear_radius_fm
            << ",\"nuclear_surface_fm\":" << cfg.nuclear_surface_fm;
    out << '}';
    return out.str();
}
