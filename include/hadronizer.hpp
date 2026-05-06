#pragma once

/**
 * @file hadronizer.hpp
 * @brief PYTHIA string-fragmentation wrapper used by eHIJING.
 */

#include "Pythia8/Pythia.h"

#include <optional>
#include <random>
#include <string>
#include <vector>

/**
 * @brief Converts medium-modified eHIJING partonic events into final hadrons.
 *
 * Hadronizer owns a dedicated PYTHIA instance configured from the
 * hadronization input file.  Each call imports the current partonic event
 * record, applies eHIJING particle-data overrides, runs string fragmentation,
 * and returns the resulting PYTHIA particles when fragmentation succeeds.
 */
class Hadronizer
{
public:
    /**
     * @brief Construct and configure the internal PYTHIA hadronizer.
     *
     * @param hadronization_config_path Path to the PYTHIA settings file used
     *        for string fragmentation and particle-data options.
     */
    explicit Hadronizer(const std::string& hadronization_config_path);

    /**
     * @brief Hadronize one medium-modified PYTHIA event record.
     *
     * @param pythiaIn Source PYTHIA generator containing the partonic event to
     *        fragment.
     * @param atomic_number Atomic number \f$ Z \f$ of the target nucleus.
     * @param mass_number Mass number \f$ A \f$ of the target nucleus.
     * @param Rx Hard-vertex x coordinate used for production vertices.
     * @param Ry Hard-vertex y coordinate used for production vertices.
     * @param Rz Hard-vertex z coordinate used for production vertices.
     *
     * @return Final PYTHIA particles if fragmentation succeeds; otherwise
     *         `std::nullopt`.
     */
    std::optional<std::vector<Pythia8::Particle>> hadronize(Pythia8::Pythia& pythiaIn,
                                             int atomic_number,
                                             int mass_number,
                                             double Rx,
                                             double Ry,
                                             double Rz);

private:
    /// Dedicated PYTHIA instance used only for hadronization.
    Pythia8::Pythia pythia;

    /// Entropy source for initializing the local random-number engine.
    std::random_device rd;

    /// Local random-number engine used by hadronization helpers.
    std::mt19937 gen;

    /// Uniform distribution on [0, 1) used by hadronization helpers.
    std::uniform_real_distribution<double> dist;
};
