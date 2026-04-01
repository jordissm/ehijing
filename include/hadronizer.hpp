#pragma once

#include "Pythia8/Pythia.h"

#include <random>
#include <vector>
#include <optional>

class Hadronizer
{
public:
    using HadronizedEvent = std::vector<Pythia8::Particle>;

    Hadronizer();

    std::optional<HadronizedEvent> hadronize(Pythia8::Pythia& pythiaIn,
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