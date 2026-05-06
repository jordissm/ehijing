#include "dis_kinematics.hpp"
#include "ehijing_constants.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
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

void require_key(bool present, const std::string& key, const std::string& path) {
    if (!present) {
        throw std::runtime_error("Missing required DIS cut `" + key + "` in " + path);
    }
}

} // namespace

DISCuts load_dis_cuts(const std::string& config_path) {
    std::ifstream input(config_path);
    if (!input) {
        throw std::runtime_error("Failed to open DIS cuts config file: " + config_path);
    }

    DISCuts cuts{};
    bool has_y_min = false;
    bool has_y_max = false;
    bool has_x_min = false;
    bool has_x_max = false;
    bool has_q2_min = false;
    bool has_w2_min = false;

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

        if (key == "yMin") {
            cuts.y_min = parsed;
            has_y_min = true;
        } else if (key == "yMax") {
            cuts.y_max = parsed;
            has_y_max = true;
        } else if (key == "xBMin") {
            cuts.bjorken_x_min = parsed;
            has_x_min = true;
        } else if (key == "xBMax") {
            cuts.bjorken_x_max = parsed;
            has_x_max = true;
        } else if (key == "Q2Min") {
            cuts.q2_min_gev2 = parsed;
            has_q2_min = true;
        } else if (key == "W2Min") {
            cuts.w2_min_gev2 = parsed;
            has_w2_min = true;
        } else {
            throw std::runtime_error(
                "Unknown DIS cut key `" + key + "` in " + config_path + ":" +
                std::to_string(line_number));
        }
    }

    require_key(has_y_min, "yMin", config_path);
    require_key(has_y_max, "yMax", config_path);
    require_key(has_x_min, "xBMin", config_path);
    require_key(has_x_max, "xBMax", config_path);
    require_key(has_q2_min, "Q2Min", config_path);
    require_key(has_w2_min, "W2Min", config_path);

    if (!(cuts.y_min < cuts.y_max)) {
        throw std::runtime_error("DIS cuts require yMin < yMax in " + config_path);
    }
    if (!(cuts.bjorken_x_min < cuts.bjorken_x_max)) {
        throw std::runtime_error("DIS cuts require xBMin < xBMax in " + config_path);
    }
    if (cuts.q2_min_gev2 < 0.0 || cuts.w2_min_gev2 < 0.0) {
        throw std::runtime_error("DIS cuts require non-negative Q2Min and W2Min in " +
                                 config_path);
    }

    return cuts;
}

bool is_valid_dis_event(const DISKinematics& kin, const DISCuts& cuts) {
    return (cuts.y_min < kin.y) &&
           (kin.y < cuts.y_max) &&
           (cuts.bjorken_x_min < kin.bjorken_x) &&
           (kin.bjorken_x < cuts.bjorken_x_max) &&
           (cuts.w2_min_gev2 < kin.W2) &&
           (cuts.q2_min_gev2 < kin.Q2);
}

bool trigger(const Pythia8::Pythia& pythia, const DISCuts& cuts) {
    return is_valid_dis_event(compute_dis_kinematics(pythia), cuts);
}
