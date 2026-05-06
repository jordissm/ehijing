#pragma once

/**
 * @file modified_ff.hpp
 * @brief Low-\f$ Q^2 \f$ medium-modified fragmentation model for eHIJING.
 */

#include "dis_kinematics.hpp"

#include "Pythia8/Pythia.h"
#include <random>
#include <string>

/**
 * @brief Monte Carlo sampler for low-\f$ Q^2 \f$ medium-induced shower changes.
 *
 * The class applies the modified fragmentation-function model after the
 * high-\f$ Q^2 \f$ shower and before hadronization.  It samples medium-induced
 * gluons and recoil partons from the eHIJING multiple-collision and nuclear
 * geometry models.
 */
class Modified_FF
{
public:
    /**
     * @brief Construct the medium-modification sampler.
     *
     * @param mode_ Medium-modification mode.  Mode 0 uses HT-style sampling;
     *        mode 1 uses generalized HT-style sampling.
     * @param atomic_number_ Atomic number \f$ Z \f$ of the target nucleus.
     * @param mass_number_ Mass number \f$ A \f$ of the target nucleus.
     * @param k_factor Overall normalization of the medium modification.
     * @param xg_n Small-\f$ x \f$ gluon-density normalization parameter.
     * @param xg_lambda Small-\f$ x \f$ gluon-density power parameter.
     * @param table_path Path used for tabulated eHIJING collision kernels.
     */
    Modified_FF(int mode_,
                int atomic_number_,
                int mass_number_,
                double k_factor,
                double xg_n,
                double xg_lambda,
                std::string table_path);

    /**
     * @brief Sample and append medium-induced partons for the current event.
     *
     * This mutates the PYTHIA event record in place and updates the hard-vertex
     * coordinates used later by hadronization/output.
     *
     * @param pythia PYTHIA generator containing the current partonic event.
     * @param kinematics Reconstructed DIS kinematics for the event.
     * @param Rx Output hard-vertex x coordinate.
     * @param Ry Output hard-vertex y coordinate.
     * @param Rz Output hard-vertex z coordinate.
     */
    void sample_ff_partons(Pythia8::Pythia& pythia,
                           const DISKinematics& kinematics,
                           double& Rx,
                           double& Ry,
                           double& Rz);

private:
    /// Selected medium-modification mode.
    int mode_, atomic_number_, mass_number_;

    /// Target proton fraction, \f$ Z/A \f$.
    double z_over_a_;

    /// eHIJING multiple-collision sampler used for medium kicks.
    EHIJING::MultipleCollision collision_sampler_;

    /// eHIJING nuclear geometry model for path lengths and thicknesses.
    EHIJING::NuclearGeometry ehijing_geometry_;

    /// Entropy source for initializing the local random-number engine.
    std::random_device random_device_;

    /// Local random-number engine for medium-modification sampling.
    std::mt19937 rng_;

    /// Uniform distribution on [0, 1) used by the sampler.
    std::uniform_real_distribution<double> uniform_dist_;
};
