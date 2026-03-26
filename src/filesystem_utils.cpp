#include "filesystem_utils.hpp"

#include <iomanip>
#include <sstream>

std::string zero_pad_int(int64_t value, int width) {
    std::ostringstream os;
    os << std::setw(width) << std::setfill('0') << value;
    return os.str();
}

std::filesystem::path shard_dir_for_event(const std::filesystem::path& base_events_dir,
                             int64_t event_id,
                             int64_t chunk_size) {
    const int64_t shard_begin = (event_id / chunk_size) * chunk_size;
    const int64_t shard_end   = shard_begin + chunk_size - 1;

    std::ostringstream name;
    name << "events_" << zero_pad_int(shard_begin)
         << "-" << zero_pad_int(shard_end);

    return base_events_dir / name.str();
}

EventPaths make_event_paths(const std::filesystem::path& base_events_dir,
                            int64_t event_id,
                            int64_t chunk_size) {

    // Set the shard directory according to the event ID and chunk size
    const std::filesystem::path shard_dir = shard_dir_for_event(base_events_dir, event_id, chunk_size);
    
    // Create the shard directory if it does not exist
    try {
        std::filesystem::create_directories(shard_dir);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "ERROR: cannot create shard directory:\n  " << shard_dir << "\n  " << e.what() << std::endl;
        return 3;
    }

    // Set the zero-padded event ID string for file naming
    const std::string event_str = zero_pad_int(event_id);

    return EventPaths{
        shard_dir,
        shard_dir / ("event_" + event_str + ".oscar"),
        shard_dir / ("event_" + event_str + ".meta.json")
    };
}