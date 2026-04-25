#include "dis_kinematics.hpp"

#include <cmath>

DISKinematics compute_dis_kinematics(const Pythia8::Pythia& pythia) {
    DISKinematics kin;

    // PYTHIA DIS event-record convention documented in dis_kinematics.hpp.
    kin.pProton = pythia.event[1].p();  // four-momentum of the incoming proton
    kin.peIn    = pythia.event[4].p();  // four-momentum of the incoming electron
    kin.peOut   = pythia.event[6].p();  // four-momentum of the outgoing electron
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
    constexpr double ymin  = 0.10;
    constexpr double ymax  = 0.85;
    constexpr double bjorken_x_min = 0.023;
    constexpr double bjorken_x_max = 0.6;
    constexpr double Q2min = 1.0;               // GeV^2
    constexpr double Wmin  = std::sqrt(10.0);   // GeV
    
    return (ymin < kin.y) && 
           (kin.y < ymax) &&
           (bjorken_x_min < kin.bjorken_x) && 
           (kin.bjorken_x < bjorken_x_max) &&
           (Wmin < kin.W) &&
           (Q2min < kin.Q2);
}

bool trigger(const Pythia8::Pythia& pythia) {
    return is_valid_dis_event(compute_dis_kinematics(pythia));
}