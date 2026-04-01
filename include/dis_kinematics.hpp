#pragma once

#include "Pythia8/Pythia.h"

struct DISKinematics {
    Pythia8::Vec4 pProton;
    Pythia8::Vec4 peIn;
    Pythia8::Vec4 peOut;
    Pythia8::Vec4 pGamma;
    double nu = 0.0;
    double Q2 = 0.0;
    double Q = 0.0;
    double W2 = 0.0;
    double W = 0.0;
    double bjorken_x = 0.0;
    double y = 0.0;
};

DISKinematics compute_dis_kinematics(const Pythia8::Pythia& pythia);
bool is_valid_dis_event(const DISKinematics& kin);
bool trigger(const Pythia8::Pythia& pythia);