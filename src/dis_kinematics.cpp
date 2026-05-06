#include "dis_kinematics.hpp"
#include "ehijing_constants.hpp"

#include <cmath>

DISKinematics compute_dis_kinematics(const Pythia8::Pythia& pythia) {
    namespace constants = ehijing::constants;

    DISKinematics kin;

    // PYTHIA DIS event-record convention documented in dis_kinematics.hpp.
    kin.pProton = pythia.event[constants::pythia::incoming_target_index].p();
    kin.peIn    = pythia.event[constants::pythia::incoming_lepton_index].p();
    kin.peOut   = pythia.event[constants::pythia::outgoing_lepton_index].p();
    kin.pGamma  = kin.peIn - kin.peOut; // four-momentum of the virtual photon/Z^0/W^±

    const double P_dot_q = kin.pProton * kin.pGamma;
    const double P_dot_k = kin.pProton * kin.peIn;

    kin.nu = P_dot_q / std::sqrt(kin.pProton * kin.pProton);
    kin.Q2 = -kin.pGamma.m2Calc();
    kin.Q  = std::sqrt(kin.Q2);
    kin.W2 = (kin.pProton + kin.pGamma).m2Calc();
    kin.W  = std::sqrt(kin.W2);
    kin.bjorken_x = kin.Q2 / (2.0 * P_dot_q);
    kin.y  = P_dot_q / P_dot_k;

    return kin;
}

bool is_valid_dis_event(const DISKinematics& kin) {
    namespace constants = ehijing::constants;
    
    return (constants::dis_cuts::y_min < kin.y) &&
           (kin.y < constants::dis_cuts::y_max) &&
           (constants::dis_cuts::bjorken_x_min < kin.bjorken_x) &&
           (kin.bjorken_x < constants::dis_cuts::bjorken_x_max) &&
           (constants::dis_cuts::w2_min_gev2 < kin.W2) &&
           (constants::dis_cuts::q2_min_gev2 < kin.Q2);
}

bool trigger(const Pythia8::Pythia& pythia) {
    return is_valid_dis_event(compute_dis_kinematics(pythia));
}
