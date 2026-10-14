#pragma once

#include "eHIJING/SeedPolicy.h"
#include <array>
#include <sstream>
#include <string>

namespace ehijing {
inline constexpr const char* seed_policy = "electra-seeds-v1";

inline std::string seed_metadata(std::uint64_t master, std::uint64_t first,
                                 std::uint64_t attempt) {
    using EHIJING::seeds::Stream;
    const std::array<const char*, 16> names = {{
        "hard_pythia", "parton_geometry", "shower_geometry", "shower_collisions",
        "shower_radiation", "modified_ff", "modified_collisions", "modified_geometry",
        "color_choice", "hadron_pythia", "hadron_helpers", "spectators",
        "hadron_parton_geometry", "hadron_shower_geometry",
        "hadron_shower_collisions", "hadron_shower_radiation"}};
    const auto context = EHIJING::seeds::chunk(master, first);
    std::ostringstream out;
    // Decimal strings retain all 64 bits in JSON readers using binary64 numbers.
    out << "{\"policy\":\"" << seed_policy << "\",\"master_seed\":\"" << master
        << "\",\"first_event_id\":" << first << ",\"attempt\":" << attempt
        << ",\"streams\":{";
    for (std::size_t i = 0; i < names.size(); ++i) {
        if (i) out << ',';
        const auto stream = static_cast<Stream>(i + 1);
        const auto value = EHIJING::seeds::derive(context, stream);
        out << '"' << names[i] << "\":";
        if (stream == Stream::HardPythia || stream == Stream::HadronPythia)
            out << EHIJING::seeds::pythia_seed(value);
        else out << '"' << value << '"';
    }
    out << "}}";
    return out.str();
}
} // namespace ehijing
