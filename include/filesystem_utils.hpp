#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

std::string zero_pad_int(int64_t value, int width = 8);

std::filesystem::path shard_dir_for_event(const std::filesystem::path& base_events_dir,
                             int64_t event_id,
                             int64_t chunk_size);

struct EventPaths {
    std::filesystem::path shardDir;
    std::filesystem::path eventPath;
    std::filesystem::path metaPath;
};

EventPaths make_event_paths(const std::filesystem::path& base_events_dir,
                            int64_t event_id,
                            int64_t chunk_size);