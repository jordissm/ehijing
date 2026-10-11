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

enum class StruckNucleon {
    Proton,
    Neutron
};

struct HadronizationResult {
    /// Momentum, energy and mass in GeV; absolute vertices in mm and coordinate
    /// times in mm/c, including breakup metadata.
    std::vector<Pythia8::Particle> particles;
    StruckNucleon struck_nucleon;
};

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
     * Rx/Ry/Rz retain the eHIJING natural-unit convention, GeV^{-1}.
     * Imported particle vertices must be displacements from the hard
     * interaction in mm (mm/c for time). The nuclear origin and reconstructed
     * displacements are converted onsite to mm and added before fragmentation.
     * The nuclear origin is added exactly once by this method.
     *
     * The source shower must have PartonVertex:setVertex disabled: that
     * subsystem assigns fm vertices in this PYTHIA version and would require
     * a separate import conversion. See the unit conventions document.
     *
     * @param pythiaIn Source PYTHIA generator containing the partonic event to
     *        fragment.
     * @param atomic_number Atomic number \f$ Z \f$ of the target nucleus.
     * @param mass_number Mass number \f$ A \f$ of the target nucleus.
     * @param Rx Hard-vertex x coordinate in GeV^{-1} used for production vertices.
     * @param Ry Hard-vertex y coordinate in GeV^{-1} used for production vertices.
     * @param Rz Hard-vertex z coordinate in GeV^{-1} used for production vertices.
     *
     * @return Hadronization result if fragmentation succeeds; otherwise
     *         `std::nullopt`.
     */
    std::optional<HadronizationResult> hadronize(Pythia8::Pythia& pythiaIn,
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
