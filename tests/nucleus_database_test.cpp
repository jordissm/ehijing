#include "cli.hpp"
#include "nucleus_database.hpp"
#include <yaml-cpp/yaml.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}
template<class F> void rejects(F action, const std::string& message) {
    try { action(); }
    catch (const std::exception& error) {
        require(std::string(error.what()).find(message) != std::string::npos, error.what());
        return;
    }
    throw std::runtime_error("Expected rejection: " + message);
}
RunConfig parse(std::vector<std::string> args) {
    std::vector<char*> argv;
    for (auto& arg : args) argv.push_back(arg.data());
    return parse_args(static_cast<int>(argv.size()), argv.data());
}
} // namespace

int main(int argc, char* argv[]) {
    // A small probe lets integration tests exercise the same parser without
    // linking the event generator or producing any simulation events.
    if (argc > 1 && std::string(argv[1]) != "--test-database") {
        try { std::cout << describe_nucleus(parse_args(argc, argv)) << '\n'; return 0; }
        catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
    }
    const auto scratch = std::filesystem::temp_directory_path() /
        ("electra-nucleus-test-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        require(argc == 3, "Expected --test-database PATH");
        const std::filesystem::path database = argv[2];
        for (const auto& item : YAML::LoadFile(database.string())["nuclei"]) {
            const auto name = item.first.as<std::string>();
            const auto target = ehijing::load_nucleus(database, name);
            require(target.name == name && target.A >= target.Z && target.nucleon_width_fm > 0., name);
        }
        const auto carbon = parse({"probe", "--nucleus", "C12", "--nucleus-database", database.string(),
                                   "--describe-nucleus"});
        const auto expected_carbon = ehijing::load_nucleus(database, "C12");
        require(carbon.atomic_number == expected_carbon.Z && carbon.mass_number == expected_carbon.A
                && carbon.nuclear_radius_fm == expected_carbon.nuclear_radius_fm,
                "C12 was not resolved");
        const auto run = parse({"probe", "--nucleus", "C12", "--nucleus-database", database.string(),
            "--number-of-events", "3", "--K", "4", "--medium-modification-mode", "1",
            "--tabulation-path", "unused", "--run-path", "unused", "--hard-process-config", "unused",
            "--hadronization-config", "unused", "--dis-cuts-config", "unused", "--seed", "12345"});
        require(!run.describe_nucleus && run.n_events == 3 && run.chunk_size == 3 && run.seed == 12345
                && run.mass_number == expected_carbon.A && run.atomic_number == expected_carbon.Z,
                "Named target was not parsed for a normal generation request");
        const auto override = parse({"probe", "--nucleus", "C12", "--nucleus-database", database.string(),
            "--nucleon-width-fm", ".7", "--nuclear-radius-fm", "2.5", "--describe-nucleus"});
        require(override.nucleon_width_fm == .7 && override.nuclear_radius_fm == 2.5,
                "Explicit override did not win");
        rejects([&] { parse({"probe", "--nucleus", "C12", "--A", "12", "--describe-nucleus"}); },
                "not both");
        rejects([&] { ehijing::load_nucleus(database, "__missing_nucleus__"); }, "available:");
        rejects([&] { ehijing::resolve_nucleus_database((scratch / "missing").string(), "probe"); },
                "not found");

        std::filesystem::create_directories(scratch);
        const auto fixture = scratch / "custom.yaml";
        const std::string valid = "nuclei:\n  d:\n    A: 2\n    Z: 1\n    distribution: hulthen\n"
            "    nucleon_width_fm: 0.61\n    hulthen_alpha_fm_inv: 0.3\n    hulthen_beta_fm_inv: 1.4\n";
        { std::ofstream out(fixture); out << valid; }
        auto custom = parse({"probe", "--nucleus", "d", "--nucleus-database", fixture.string(),
                             "--describe-nucleus"});
        require(custom.nucleon_width_fm == .61 && custom.hulthen_alpha_fm_inv == .3
                && custom.hulthen_beta_fm_inv == 1.4, "Hulthen values were ignored");
        rejects([&] { parse({"probe", "--nucleus", "d", "--nucleus-database", fixture.string(),
            "--nuclear-radius-fm", "2", "--describe-nucleus"}); }, "do not apply");
        for (const auto& replacement : std::vector<std::pair<std::string, std::string>>{
                {"A: 2", "A: 3"}, {"0.61", ".nan"}, {"0.3", "1.5"}, {"Z: 1", "Z: 1\n    Z: 1"},
                {"distribution: hulthen", "distribution: unknown"}, {"0.61", "-0.5"},
                {"nucleon_width_fm", "misspelled_width"}}) {
            auto invalid = valid;
            invalid.replace(invalid.find(replacement.first), replacement.first.size(), replacement.second);
            { std::ofstream out(fixture); out << invalid; }
            rejects([&] { ehijing::load_nucleus(fixture, "d"); }, "Nucleus 'd'");
        }
        { std::ofstream out(fixture); out << "nuclei: ["; }
        rejects([&] { ehijing::load_nucleus(fixture, "d"); }, "Nucleus 'd'");
        std::filesystem::remove_all(scratch);
        std::cout << "Nucleus lookup, validation, overrides and Hulthen configuration passed\n";
    } catch (const std::exception& error) {
        std::filesystem::remove_all(scratch);
        std::cerr << error.what() << '\n'; return 1;
    }
}
