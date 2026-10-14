#pragma once
#include <filesystem>
#include <string>

namespace ehijing {
// SHA-256 of file bytes. No path, timestamp, locale or RNG enters the digest.
std::string file_sha256(const std::filesystem::path& path);
}
