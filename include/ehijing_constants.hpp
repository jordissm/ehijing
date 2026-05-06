#pragma once

/**
 * @file ehijing_constants.hpp
 * @brief Central compile-time constants shared by the eHIJING driver.
 */

#include <cstdint>

/// Compile-time constants shared across the eHIJING generator components.
namespace ehijing::constants {

/// Mathematical constants.
namespace math {
/// Ratio of a circle circumference to its diameter.
inline constexpr double pi = 3.141592653589793238462643383279502884;

/// Full azimuthal angle in radians.
inline constexpr double two_pi = 2.0 * pi;
} // namespace math

/// Unit-conversion factors.
namespace units {
/// \f$ \hbar c \f$ in GeV fm.
inline constexpr double hbarc_gev_fm = 0.19732698;
} // namespace units

/// PDG Monte Carlo particle IDs used by the eHIJING workflow.
///
/// Particle constants in this namespace intentionally use the `_id` suffix to
/// distinguish PDG identifiers from masses, statuses, and array indices.
namespace pdg {
inline constexpr int down_id = 1;
inline constexpr int up_id = 2;
inline constexpr int strange_id = 3;
inline constexpr int gluon_id = 21;

inline constexpr int pi0_id = 111;
inline constexpr int pion_charged_id = 211;
inline constexpr int kaon_charged_id = 321;
inline constexpr int k0_id = 311;
inline constexpr int k0_long_id = 130;
inline constexpr int k0_short_id = 310;

inline constexpr int neutron_id = 2112;
inline constexpr int proton_id = 2212;

inline constexpr int dd_spin1_diquark_id = 1103;
inline constexpr int ud_spin0_diquark_id = 2101;
inline constexpr int ud_spin1_diquark_id = 2103;
inline constexpr int uu_spin1_diquark_id = 2203;

/// Lower exclusive absolute-PDG boundary used to identify diquarks.
inline constexpr int diquark_id_min_exclusive = 1000;

/// Upper exclusive absolute-PDG boundary used to identify diquarks.
inline constexpr int diquark_id_max_exclusive = 3000;
} // namespace pdg

/// Mass values used when synchronizing PYTHIA, eHIJING, and output records.
namespace mass {
/// Effective charged/neutral pion mass in GeV for particle-data overrides.
inline constexpr double pion_gev = 0.138;

/// Effective kaon mass in GeV for particle-data overrides.
inline constexpr double kaon_gev = 0.494;

/// Proton mass in GeV used for particle-data overrides.
inline constexpr double proton_gev = 0.938;

/// Neutron mass in GeV used for particle-data overrides.
inline constexpr double neutron_gev = 0.938;

/// Constituent light-quark mass in GeV used by remnant construction.
inline constexpr double constituent_light_quark_gev = 0.33;

/// Remnant diquark mass in GeV used by recoil construction.
inline constexpr double remnant_diquark_gev = 0.57933;
} // namespace mass

/// Nuclear-geometry constants.
namespace nuclear {
/// Sharp-radius coefficient in \f$ R_A = r_0 A^{1/3} \f$, in fm.
inline constexpr double radius_coefficient_fm = 1.2;
} // namespace nuclear

/// Small numerical cutoffs used to guard rotations and kinematic inversions.
namespace numeric {
/// Tolerance for deciding whether a momentum rotation is needed.
inline constexpr double rotation_epsilon = 1.0e-6;

/// Minimum accumulated transverse momentum kick squared.
inline constexpr double min_momentum_kick2 = 1.0e-9;

/// Minimum allowed \f$ \sin^2\theta \f$ in frame rotations.
inline constexpr double min_sin_theta2 = 1.0e-9;

/// Minimum transverse momentum scale for daughter construction.
inline constexpr double min_kt_daughter = 1.0e-10;

/// Small positive formation-time offset used for output vertices, in fm.
inline constexpr double formation_epsilon_fm = 1.0e-4;
} // namespace numeric

/// PYTHIA event-record conventions and integer settings used by eHIJING.
namespace pythia {
/// Maximum seed accepted by PYTHIA's random-number settings.
inline constexpr std::uint64_t max_seed = 900000000ULL;

/// Incoming target index in the PYTHIA DIS event record.
inline constexpr int incoming_target_index = 1;

/// Incoming lepton index in the PYTHIA DIS event record.
inline constexpr int incoming_lepton_index = 4;

/// Outgoing hard-parton index in the PYTHIA DIS event record.
inline constexpr int hard_parton_index = 5;

/// Outgoing scattered-lepton index in the PYTHIA DIS event record.
inline constexpr int outgoing_lepton_index = 6;

/// PYTHIA status code for outgoing hard-process particles.
inline constexpr int hard_process_outgoing_status = 23;

/// PYTHIA status code for beam remnants.
inline constexpr int beam_remnant_status = 63;

/// Custom status code for medium recoil partons.
inline constexpr int medium_recoil_status = 201;

/// PDG base used by PYTHIA ion IDs.
inline constexpr int nucleus_pdg_base = 100000000;

/// Multiplicative factor for the nuclear charge in a PYTHIA ion ID.
inline constexpr int nucleus_pdg_z_factor = 10000;

/// Multiplicative factor for the mass number in a PYTHIA ion ID.
inline constexpr int nucleus_pdg_a_factor = 10;
} // namespace pythia

/// Constants specific to the modified-fragmentation-function stage.
namespace modified_ff {
/// Minimum parton energy retained in medium-modification sampling, in GeV.
inline constexpr double min_parton_energy_gev = 0.2;
} // namespace modified_ff

} // namespace ehijing::constants
