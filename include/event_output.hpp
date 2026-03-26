#pragma once

#include "Pythia8/Pythia.h"
#include "dis_kinematics.hpp"

#include <cstdint>
#include <iosfwd>
#include <vector>

void write_event_headers(std::ostream& out);

void write_event_output(
    int32_t event_number,
    int Z,
    int A,
    const DISKinematics& kin,
    const std::vector<Pythia8::Particle>& particles,
    std::ostream& event_out,
    std::ostream& meta_out);