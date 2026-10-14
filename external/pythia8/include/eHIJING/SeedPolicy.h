#ifndef EHIJING_SEED_POLICY_H
#define EHIJING_SEED_POLICY_H

#include <cstdint>
#include <random>

namespace EHIJING {
namespace seeds {

// Frozen identifiers and unsigned arithmetic define electra-seeds-v1.
// Changing either requires a new policy version and new reference vectors.
enum class Stream : std::uint64_t {
  HardPythia = 1, PartonGeometry = 2, ShowerGeometry = 3,
  ShowerCollisions = 4, ShowerRadiation = 5, ModifiedFF = 6,
  ModifiedCollisions = 7, ModifiedGeometry = 8, ColorChoice = 9,
  HadronPythia = 10, HadronHelpers = 11, Spectators = 12,
  HadronPartonGeometry = 13, HadronShowerGeometry = 14,
  HadronShowerCollisions = 15, HadronShowerRadiation = 16
};

inline std::uint64_t mix(std::uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

inline std::uint64_t chunk(std::uint64_t master, std::uint64_t first) {
  return mix(master ^ mix(first) ^ UINT64_C(0x454c454354524131));
}

inline std::uint64_t derive(std::uint64_t context, Stream stream) {
  return mix(context ^ mix(static_cast<std::uint64_t>(stream)));
}

inline int pythia_seed(std::uint64_t value) {
  return static_cast<int>(1 + value % UINT64_C(900000000));
}

inline void seed(std::mt19937& engine, std::uint64_t value) {
  std::seed_seq words{static_cast<std::uint32_t>(value),
                      static_cast<std::uint32_t>(value >> 32)};
  engine.seed(words);
}

inline std::mt19937 engine(std::uint64_t value) {
  std::mt19937 result;
  seed(result, value);
  return result;
}

} // namespace seeds

// Linked symbol: an application built for this policy cannot silently use an
// older medium library whose constructors still obtain independent entropy.
int seed_policy_version();
} // namespace EHIJING
#endif
