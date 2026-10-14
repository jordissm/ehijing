#include "random_seed.hpp"
#include "file_digest.hpp"
#include "cli.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    auto path = std::filesystem::temp_directory_path() /
        ("electra-digest-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        using namespace EHIJING::seeds;
        require(EHIJING::seed_policy_version() == 1, "Old medium library");
        require(chunk(0, 0) == UINT64_C(13059828666866575528), "Zero-seed context changed");
        const auto context = chunk(12345, 0);
        require(context == UINT64_C(5504340652873255536), "Chunk vector changed");
        require(derive(context, Stream::ModifiedFF) == UINT64_C(8169500210558457899), "Stream vector changed");
        require(pythia_seed(derive(context, Stream::HardPythia)) == 91652349, "Hard PYTHIA vector changed");
        require(pythia_seed(derive(context, Stream::HadronPythia)) == 742193345, "Hadron PYTHIA vector changed");
        require(chunk(UINT64_MAX, 123) == UINT64_C(9650306775149182456), "Wide-seed context changed");
        require(pythia_seed(0) == 1 && pythia_seed(899999999) == 900000000
            && pythia_seed(900000000) == 1, "PYTHIA boundary conversion failed");
        auto first = engine(17), replay = engine(17), wide = engine(UINT64_C(1) << 32 | 17);
        bool differs = false;
        for (int i = 0; i < 100; ++i) {
            const auto draw = first();
            require(draw == replay(), "Engine replay failed");
            differs |= draw != wide();
        }
        require(differs, "High master-seed bits lost");
        std::vector<std::string> args = {"ehijing", "--number-of-events", "1", "--Z", "1",
            "--A", "2", "--K", "4", "--medium-modification-mode", "0",
            "--tabulation-path", "unused", "--run-path", "unused",
            "--hard-process-config", "unused", "--hadronization-config", "unused",
            "--dis-cuts-config", "unused", "--seed", "18446744073709551615"};
        std::vector<char*> argv;
        for (auto& arg : args) argv.push_back(arg.data());
        require(parse_args(argv.size(), argv.data()).seed == UINT64_MAX, "CLI reduced master seed");
        std::ofstream(path, std::ios::binary).close();
        require(ehijing::file_sha256(path) == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
            "Empty SHA-256 vector failed");
        { std::ofstream file(path, std::ios::binary); file << "abc"; }
        require(ehijing::file_sha256(path) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
            "Short SHA-256 vector failed");
        { std::ofstream file(path, std::ios::binary); for (int i=0;i<1000;++i) file << 'a'; }
        require(ehijing::file_sha256(path) == "41edece42d63e8d9bf515a9ba6932e1c20cbc9f5a5d134645adb5db1b9737ea3",
            "Multiblock SHA-256 vector failed");
        std::filesystem::remove(path);
        std::cout << "Seed policy, parser and digest checks passed\n";
    } catch (const std::exception& error) {
        std::filesystem::remove(path);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
