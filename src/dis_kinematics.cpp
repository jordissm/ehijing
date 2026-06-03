#include "dis_kinematics.hpp"
#include "ehijing_constants.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>

DISKinematics compute_dis_kinematics(const Pythia8::Pythia& pythia) {
    namespace constants = ehijing::constants;

    DISKinematics kin;

    // PYTHIA DIS event-record convention documented in dis_kinematics.hpp.
    kin.pProton = pythia.event[constants::pythia::incoming_target_index].p();
    kin.peIn    = pythia.event[constants::pythia::incoming_lepton_index].p();
    kin.peOut   = pythia.event[constants::pythia::outgoing_lepton_index].p();
    kin.pGamma  = kin.peIn - kin.peOut; // four-momentum of the virtual photon/Z^0/W^±

    const double P_dot_q = kin.pProton * kin.pGamma;
    const double P_dot_k = kin.pProton * kin.peIn;

    kin.nu = P_dot_q / std::sqrt(kin.pProton * kin.pProton);
    kin.Q2 = -kin.pGamma.m2Calc();
    kin.Q  = std::sqrt(kin.Q2);
    kin.W2 = (kin.pProton + kin.pGamma).m2Calc();
    kin.W  = std::sqrt(kin.W2);
    kin.bjorken_x = kin.Q2 / (2.0 * P_dot_q);
    kin.y  = P_dot_q / P_dot_k;

    return kin;
}

namespace {

std::string trim(std::string text) {
    auto is_not_space = [](unsigned char c) { return !std::isspace(c); };

    text.erase(text.begin(), std::find_if(text.begin(), text.end(), is_not_space));
    text.erase(std::find_if(text.rbegin(), text.rend(), is_not_space).base(), text.end());
    return text;
}

double parse_double(const std::string& value,
                    const std::string& key,
                    const std::string& config_path,
                    int line_number) {
    std::size_t parsed_length = 0;
    try {
        const double parsed = std::stod(value, &parsed_length);
        if (trim(value.substr(parsed_length)).empty()) {
            return parsed;
        }
    } catch (const std::exception&) {
    }

    throw std::runtime_error(
        "Invalid value for " + key + " in " + config_path + ":" +
        std::to_string(line_number) + ": " + value);
}

void validate_range(const std::optional<double>& min_value,
                    const std::optional<double>& max_value,
                    const std::string& min_key,
                    const std::string& max_key,
                    const std::string& path) {
    if (min_value && max_value && !(*min_value < *max_value)) {
        throw std::runtime_error(
            "DIS cuts require " + min_key + " < " + max_key + " in " + path);
    }
}

void validate_non_negative(const std::optional<double>& value,
                           const std::string& key,
                           const std::string& path) {
    if (value && *value < 0.0) {
        throw std::runtime_error("DIS cuts require non-negative " + key + " in " + path);
    }
}

bool passes_min(double value, const std::optional<double>& min_value) {
    return !min_value || (*min_value < value);
}

bool passes_max(double value, const std::optional<double>& max_value) {
    return !max_value || (value < *max_value);
}

bool passes_range(double value,
                  const std::optional<double>& min_value,
                  const std::optional<double>& max_value) {
    return passes_min(value, min_value) && passes_max(value, max_value);
}

bool assign_cut(const std::string& key, double parsed, DISCuts& cuts) {
    if (key == "yMin") {
        cuts.y_min = parsed;
    } else if (key == "yMax") {
        cuts.y_max = parsed;
    } else if (key == "xBMin") {
        cuts.bjorken_x_min = parsed;
    } else if (key == "xBMax") {
        cuts.bjorken_x_max = parsed;
    } else if (key == "nuMin") {
        cuts.nu_min_gev = parsed;
    } else if (key == "nuMax") {
        cuts.nu_max_gev = parsed;
    } else if (key == "Q2Min") {
        cuts.q2_min_gev2 = parsed;
    } else if (key == "Q2Max") {
        cuts.q2_max_gev2 = parsed;
    } else if (key == "W2Min") {
        cuts.w2_min_gev2 = parsed;
    } else if (key == "W2Max") {
        cuts.w2_max_gev2 = parsed;
    } else {
        return false;
    }
    return true;
}

} // namespace

DISCuts load_dis_cuts(const std::string& config_path) {
    std::ifstream input(config_path);
    if (!input) {
        throw std::runtime_error("Failed to open DIS cuts config file: " + config_path);
    }

    DISCuts cuts{};

    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;

        const std::size_t comment_position = line.find('#');
        if (comment_position != std::string::npos) {
            line.erase(comment_position);
        }

        line = trim(line);
        if (line.empty()) {
            continue;
        }

        const std::size_t separator_position = line.find('=');
        if (separator_position == std::string::npos) {
            throw std::runtime_error(
                "Expected `key = value` in " + config_path + ":" +
                std::to_string(line_number));
        }

        const std::string key = trim(line.substr(0, separator_position));
        const std::string value = trim(line.substr(separator_position + 1));
        const double parsed = parse_double(value, key, config_path, line_number);

        if (!assign_cut(key, parsed, cuts)) {
            throw std::runtime_error(
                "Unknown DIS cut key `" + key + "` in " + config_path + ":" +
                std::to_string(line_number));
        }
    }

    validate_range(cuts.y_min, cuts.y_max, "yMin", "yMax", config_path);
    validate_range(cuts.bjorken_x_min, cuts.bjorken_x_max, "xBMin", "xBMax", config_path);
    validate_range(cuts.nu_min_gev, cuts.nu_max_gev, "nuMin", "nuMax", config_path);
    validate_range(cuts.q2_min_gev2, cuts.q2_max_gev2, "Q2Min", "Q2Max", config_path);
    validate_range(cuts.w2_min_gev2, cuts.w2_max_gev2, "W2Min", "W2Max", config_path);

    validate_non_negative(cuts.q2_min_gev2, "Q2Min", config_path);
    validate_non_negative(cuts.q2_max_gev2, "Q2Max", config_path);
    validate_non_negative(cuts.w2_min_gev2, "W2Min", config_path);
    validate_non_negative(cuts.w2_max_gev2, "W2Max", config_path);

    return cuts;
}

bool is_valid_dis_event(const DISKinematics& kin, const DISCuts& cuts) {
    return passes_range(kin.y, cuts.y_min, cuts.y_max) &&
           passes_range(kin.bjorken_x, cuts.bjorken_x_min, cuts.bjorken_x_max) &&
           passes_range(kin.nu, cuts.nu_min_gev, cuts.nu_max_gev) &&
           passes_range(kin.Q2, cuts.q2_min_gev2, cuts.q2_max_gev2) &&
           passes_range(kin.W2, cuts.w2_min_gev2, cuts.w2_max_gev2);
}

bool trigger(const Pythia8::Pythia& pythia, const DISCuts& cuts) {
    return is_valid_dis_event(compute_dis_kinematics(pythia), cuts);
}
