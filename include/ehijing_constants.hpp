#pragma once

#include <cstdint>

namespace ehijing::constants {

namespace math {
inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double two_pi = 2.0 * pi;
} // namespace math

namespace units {
inline constexpr double hbarc_gev_fm = 0.19732698;
} // namespace units

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
inline constexpr int diquark_id_min_exclusive = 1000;
inline constexpr int diquark_id_max_exclusive = 3000;
} // namespace pdg

namespace mass {
inline constexpr double pion_gev = 0.138;
inline constexpr double kaon_gev = 0.494;
inline constexpr double output_nucleon_gev = 0.938;
inline constexpr double proton_gev = 0.93847;
inline constexpr double neutron_gev = 0.93957;
inline constexpr double constituent_light_quark_gev = 0.33;
inline constexpr double remnant_diquark_gev = 0.57933;
} // namespace mass

namespace nuclear {
inline constexpr double radius_coefficient_fm = 1.2;
} // namespace nuclear

namespace numeric {
inline constexpr double rotation_epsilon = 1.0e-6;
inline constexpr double min_momentum_kick2 = 1.0e-9;
inline constexpr double min_sin_theta2 = 1.0e-9;
inline constexpr double min_kt_daughter = 1.0e-10;
inline constexpr double formation_epsilon_fm = 1.0e-4;
} // namespace numeric

namespace pythia {
inline constexpr std::uint64_t max_seed = 900000000ULL;

inline constexpr int incoming_target_index = 1;
inline constexpr int incoming_lepton_index = 4;
inline constexpr int hard_parton_index = 5;
inline constexpr int outgoing_lepton_index = 6;

inline constexpr int hard_process_outgoing_status = 23;
inline constexpr int beam_remnant_status = 63;
inline constexpr int medium_recoil_status = 201;

inline constexpr int nucleus_pdg_base = 100000000;
inline constexpr int nucleus_pdg_z_factor = 10000;
inline constexpr int nucleus_pdg_a_factor = 10;
} // namespace pythia

namespace modified_ff {
inline constexpr double min_parton_energy_gev = 0.2;
} // namespace modified_ff

} // namespace ehijing::constants
