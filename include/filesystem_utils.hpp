#pragma once

/**
 * @file filesystem_utils.hpp
 * @brief Helpers for deterministic, sharded eHIJING event output paths.
 */

#include <cstdint>
#include <filesystem>
#include <string>

/**
 * @brief Convert an integer to a zero-padded decimal string.
 *
 * @param value Integer value to format.
 * @param width Minimum output width. Shorter values are padded with leading
 *        zeroes.
 *
 * @return Zero-padded decimal representation of `value`.
 */
std::string zero_pad_int(int64_t value, int width = 8);

/**
 * @brief Compute the shard directory containing a given event.
 *
 * Events are grouped into directories named
 * `events_<first-event>-<last-event>`, where both event IDs are zero-padded.
 *
 * @param base_events_dir Parent directory for all event shards.
 * @param event_id Event ID whose shard should be computed.
 * @param chunk_size Number of events per shard directory.
 *
 * @return Filesystem path to the shard directory for `event_id`.
 */
std::filesystem::path shard_dir_for_event(const std::filesystem::path& base_events_dir,
                             int64_t event_id,
                             int64_t chunk_size);

/**
 * @brief Concrete output paths for one generated event.
 */
struct EventPaths {
    /// Directory containing this event's output files.
    std::filesystem::path shardDir;

    /// OSCAR2013 particle-list output path.
    std::filesystem::path eventPath;

    /// Per-event JSON metadata output path.
    std::filesystem::path metaPath;
};

/**
 * @brief Create and return all output paths for one event.
 *
 * The shard directory is created if it does not already exist.  On filesystem
 * errors the function reports the failure and exits with status code 2.
 *
 * @param base_events_dir Parent directory for all event shards.
 * @param event_id Event ID used for sharding and file naming.
 * @param chunk_size Number of events per shard directory.
 *
 * @return Paths for the shard directory, OSCAR output, and metadata output.
 */
EventPaths make_event_paths(const std::filesystem::path& base_events_dir,
                            int64_t event_id,
                            int64_t chunk_size);
