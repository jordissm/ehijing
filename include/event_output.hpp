#pragma once

/**
 * @file event_output.hpp
 * @brief Writers for eHIJING OSCAR2013 event records and metadata.
 */

#include "Pythia8/Pythia.h"
#include "dis_kinematics.hpp"

#include <cstdint>
#include <iosfwd>
#include <vector>

/**
 * @brief Write the OSCAR2013 particle-list header lines.
 *
 * @param out Stream that receives the event-file header.
 */
void write_event_headers(std::ostream& out);

/**
 * @brief Write one generated event and its metadata.
 *
 * The event stream receives an OSCAR2013-style particle list containing final
 * hadrons and spectator nucleons.  The metadata stream receives a compact JSON
 * record with the event ID, target, and reconstructed DIS kinematics.
 *
 * @param event_number Global event ID.
 * @param atomic_number Atomic number \f$ Z \f$ of the target nucleus.
 * @param mass_number Mass number \f$ A \f$ of the target nucleus.
 * @param kin Reconstructed DIS kinematics for the event.
 * @param particles PYTHIA particles returned by the hadronizer.
 * @param event_out Stream that receives the OSCAR2013 particle list.
 * @param meta_out Stream that receives the JSON metadata record.
 */
void write_event_output(
    int32_t event_number,
    int atomic_number,
    int mass_number,
    const DISKinematics& kin,
    const std::vector<Pythia8::Particle>& particles,
    std::ostream& event_out,
    std::ostream& meta_out);
