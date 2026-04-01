#include "dis_kinematics.hpp"

#include <cmath>

DISKinematics compute_dis_kinematics(const Pythia8::Pythia& pythia) {
    DISKinematics k;
    k.pProton = pythia.event[1].p(); // four-momentum of the incoming proton
    k.peIn    = pythia.event[4].p(); // four-momentum of the incoming electron
    k.peOut   = pythia.event[6].p(); // four-momentum of the outgoing electron
    k.pGamma  = k.peIn - k.peOut;    // four-momentum of the virtual photon/Z^0/W^±

    k.nu = (k.pProton * k.pGamma) / std::sqrt(k.pProton * k.pProton);
    k.Q2 = -k.pGamma.m2Calc();
    k.Q  = std::sqrt(k.Q2);
    k.W2 = (k.pProton + k.pGamma).m2Calc();
    k.W  = std::sqrt(k.W2);
    k.bjorken_x = k.Q2 / (2.0 * (k.pProton * k.pGamma));
    k.y  = (k.pProton * k.pGamma) / (k.pProton * k.peIn);

    return k;
}

bool is_valid_dis_event(const DISKinematics& kin) {
    constexpr double ymin  = 0.10;
    constexpr double ymax  = 0.85;
    constexpr double bjorken_x_min = 0.023;
    constexpr double bjorken_x_max = 0.6;
    constexpr double Q2min = 1.0; // GeV^2
    constexpr double Wmin  = std::sqrt(10.0); // GeV
    
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