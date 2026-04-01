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
                int atomic_number_,
                int mass_number_,
                double k_factor,
                double xg_n,
                double xg_lambda,
                std::string table_path);

    void sample_ff_partons(Pythia& pythia, double& Rx, double& Ry, double& Rz);

private:
    int mode_, atomic_number_, mass_number_;
    double z_over_a_;
    EHIJING::MultipleCollision collision_sampler_;
    EHIJING::NuclearGeometry ehijing_geometry_;
    std::random_device random_device_;
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_dist_;
};