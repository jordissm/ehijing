#pragma once

#include "Pythia8/Pythia.h"
#include <random>
#include <vector>

using namespace Pythia8;

class Hadronizer
{
public:
    Hadronizer();

    std::vector<Particle> hadronize(Pythia& pythiaIn,
                                    int atomic_number,
                                    int mass_number,
                                    double Rx,
                                    double Ry,
                                    double Rz);

private:
    Pythia pythia;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;
};