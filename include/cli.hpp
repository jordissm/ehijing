#pragma once

#include <cstdint>
#include <string>

struct RunConfig {
    int n_events;
    int atomic_number;
    int mass_number;
    int mode;
    double k_factor;
    std::string table_path;
    std::string runDir;
    std::string configFile;
    int64_t firstEventId;
    int64_t chunk_size;
    uint32_t seed;
};

void usage(const char* prog);
RunConfig parse_args(int argc, char* argv[]);