#pragma once

#include "Pythia8/Pythia.h"
#include <random>
#include <string>

using namespace Pythia8;

// low-Q2 medium correction (a Monte Carlo version of the the modified FF model)
class Modified_FF
{
public:
    Modified_FF(int mode_,
                int Z_,
                int A_,
                double K,
                double n,
                double lambda,
                std::string TablePath);

    void sample_FF_partons(Event& event, double& Rx, double& Ry, double& Rz);

private:
    int mode, Z, A;
    double ZoverA;
    EHIJING::MultipleCollision Coll;
    EHIJING::NuclearGeometry eHIJING_Geometry;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;
};