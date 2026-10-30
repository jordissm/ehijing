#pragma once

#include <filesystem>
#include <string>

namespace ehijing {

struct NucleusParameters {
    std::string name, distribution;
    int A = 0, Z = 0;
    double nucleon_width_fm = 0.;
    double nuclear_radius_fm = 0., nuclear_surface_fm = 0.;
    double hulthen_alpha_fm_inv = .228, hulthen_beta_fm_inv = 1.18;
};

std::filesystem::path resolve_nucleus_database(
    const std::string& requested, const std::string& executable);
NucleusParameters load_nucleus(const std::filesystem::path& path, const std::string& name);

} // namespace ehijing
