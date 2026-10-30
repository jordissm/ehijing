#include "nucleus_database.hpp"

#include <yaml-cpp/yaml.h>
#include <cmath>
#include <cstdlib>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifndef EHIJING_SOURCE_NUCLEUS_DATABASE
#define EHIJING_SOURCE_NUCLEUS_DATABASE ""
#endif
#ifndef EHIJING_INSTALLED_NUCLEUS_DATABASE
#define EHIJING_INSTALLED_NUCLEUS_DATABASE ""
#endif

namespace ehijing {
namespace {

void unique_keys(const YAML::Node& node, const std::string& context) {
    if (!node.IsMap()) throw std::runtime_error(context + " must be a mapping");
    std::set<std::string> keys;
    for (const auto& item : node) {
        const auto key = item.first.as<std::string>();
        if (!keys.insert(key).second)
            throw std::runtime_error("Duplicate key '" + key + "' in " + context);
    }
}

template<class T> T required(const YAML::Node& node, const std::string& key) {
    if (!node[key] || !node[key].IsScalar())
        throw std::runtime_error("Missing scalar field '" + key + "'");
    return node[key].as<T>();
}

double positive(const YAML::Node& node, const std::string& key) {
    const double value = required<double>(node, key);
    if (!std::isfinite(value) || value <= 0.)
        throw std::runtime_error("Field '" + key + "' must be finite and positive");
    return value;
}

} // namespace

std::filesystem::path resolve_nucleus_database(
    const std::string& requested, const std::string& executable) {
    namespace fs = std::filesystem;
    auto checked = [](const fs::path& path) {
        if (!fs::is_regular_file(path))
            throw std::runtime_error("Nucleus database not found: " + path.string());
        return fs::canonical(path);
    };
    if (!requested.empty()) return checked(requested);
    if (const char* value = std::getenv("ELECTRA_NUCLEUS_DATABASE"))
        if (*value) return checked(value);

    std::vector<fs::path> candidates = {"input/ehijing/woods_saxon.yaml"};
    fs::path binary(executable);
    if (!binary.has_parent_path()) {
        if (const char* search = std::getenv("PATH")) {
            std::istringstream paths(search);
            std::string directory;
            while (std::getline(paths, directory, ':')) {
                const fs::path candidate = fs::path(directory) / binary;
                if (fs::is_regular_file(candidate)) { binary = candidate; break; }
            }
        }
    }
    if (fs::is_regular_file(binary)) {
        const auto directory = fs::canonical(binary).parent_path();
        candidates.push_back(directory / "../share/electra/ehijing/woods_saxon.yaml");
    }
    candidates.emplace_back(EHIJING_SOURCE_NUCLEUS_DATABASE);
    candidates.emplace_back(EHIJING_INSTALLED_NUCLEUS_DATABASE);
    for (const auto& candidate : candidates)
        if (!candidate.empty() && fs::is_regular_file(candidate)) return checked(candidate);
    throw std::runtime_error("Cannot locate woods_saxon.yaml; supply --nucleus-database PATH "
                             "or set ELECTRA_NUCLEUS_DATABASE");
}

NucleusParameters load_nucleus(const std::filesystem::path& path, const std::string& name) {
    try {
        const auto root = YAML::LoadFile(path.string());
        unique_keys(root, "database");
        const auto nuclei = root["nuclei"];
        unique_keys(nuclei, "nuclei");
        const auto entry = nuclei[name];
        if (!entry) {
            std::string available;
            for (const auto& item : nuclei)
                available += (available.empty() ? "" : ", ") + item.first.as<std::string>();
            throw std::runtime_error("Unknown nucleus '" + name + "'; available: " + available);
        }
        unique_keys(entry, "nucleus '" + name + "'");
        NucleusParameters result;
        result.name = name;
        result.A = required<int>(entry, "A");
        result.Z = required<int>(entry, "Z");
        if (result.A < 1 || result.Z < 0 || result.Z > result.A)
            throw std::runtime_error("Require A >= 1 and 0 <= Z <= A");
        result.distribution = required<std::string>(entry, "distribution");
        result.nucleon_width_fm = positive(entry, "nucleon_width_fm");
        std::set<std::string> fields = {"A", "Z", "distribution", "nucleon_width_fm"};
        if (result.distribution == "woods_saxon") {
            if (result.A <= 2) throw std::runtime_error("Woods-Saxon requires A > 2");
            result.nuclear_radius_fm = positive(entry, "nuclear_radius_fm");
            result.nuclear_surface_fm = positive(entry, "nuclear_surface_fm");
            fields.insert("nuclear_radius_fm"); fields.insert("nuclear_surface_fm");
        } else if (result.distribution == "hulthen") {
            if (result.A != 2 || result.Z != 1)
                throw std::runtime_error("Hulthen requires A=2, Z=1");
            result.hulthen_alpha_fm_inv = positive(entry, "hulthen_alpha_fm_inv");
            result.hulthen_beta_fm_inv = positive(entry, "hulthen_beta_fm_inv");
            if (result.hulthen_alpha_fm_inv >= result.hulthen_beta_fm_inv)
                throw std::runtime_error("Hulthen requires alpha < beta");
            fields.insert("hulthen_alpha_fm_inv"); fields.insert("hulthen_beta_fm_inv");
        } else throw std::runtime_error("Unsupported distribution '" + result.distribution + "'");
        for (const auto& item : entry) {
            const auto key = item.first.as<std::string>();
            if (!fields.count(key)) throw std::runtime_error("Unexpected field '" + key + "'");
        }
        return result;
    } catch (const std::exception& error) {
        throw std::runtime_error("Nucleus '" + name + "' in " + path.string() + ": " + error.what());
    }
}

} // namespace ehijing
