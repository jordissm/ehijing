#pragma once

#include <cstdint>
#include <string>

struct RunConfig {
    int nEvents;
    int Z;
    int A;
    int mode;
    double K;
    std::string tableDir;
    std::string runDir;
    std::string configFile;
    int64_t firstEventId;
    int64_t chunkSize;
    uint32_t seed;
};

void usage(const char* prog);
RunConfig parse_args(int argc, char* argv[]);