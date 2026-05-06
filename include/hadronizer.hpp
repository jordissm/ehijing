#pragma once

#include "Pythia8/Pythia.h"

#include <random>
#include <string>
#include <vector>
#include <optional>

class Hadronizer
{
public:
    explicit Hadronizer(const std::string& hadronization_config_path);

    std::optional<std::vector<Pythia8::Particle>> hadronize(Pythia8::Pythia& pythiaIn,
                                             int atomic_number,
                                             int mass_number,
                                             double Rx,
                                             double Ry,
                                             double Rz);

private:
    Pythia8::Pythia pythia;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;
};
