#pragma once

/**
 * @file event_output.hpp
 * @brief Writers for eHIJING OSCAR2013 event records and metadata.
 */

#include "Pythia8/Pythia.h"
#include "dis_kinematics.hpp"
#include "hadronizer.hpp"

#include <cstdint>
#include <iosfwd>
#include <vector>

class SpacetimeHistory;

/** Spectator counts required by neutral-current DIS baryon/charge balance.
 * Quantum numbers and required counts are retained in thirds to avoid
 * truncating fractional quantum numbers. Legacy mode keeps the A-1 fallback;
 * explicit-target mode rejects any mismatch with untouched identities.
 */
struct SpectatorBalance {
    int64_t final_baryon_number_times3 = 0;
    int64_t final_charge_times3 = 0;
    int64_t required_protons_times3 = 0;
    int64_t required_neutrons_times3 = 0;
    int64_t protons = 0;
    int64_t neutrons = 0;
    bool admissible = false;
    bool baryon_number_conserved = false;
    bool electric_charge_conserved = false;
};

/** Count only final hadronizer products, including antibaryons, resonances,
 * charged mesons and charged decay products. The scattered DIS lepton is not
 * part of this collection. Admissible required counts are nonnegative integers.
 */
SpectatorBalance calculate_spectator_balance(
    int atomic_number, int mass_number, StruckNucleon struck_nucleon,
    const std::vector<Pythia8::Particle>& particles);

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
 * Spectators satisfy B_final + N_p + N_n = A and Q_final + N_p = Z when
 * the required counts are admissible. Otherwise legacy spectator counts are
 * retained and metadata/history explicitly report the unsatisfied balance in
 * legacy mode. With an explicit target, untouched identities must match the
 * balance; a mismatch is an error and must be rejected before invoking output.
 *
 * Hadron positions are freestreamed from the saved vertex to t = 0. With
 * HadronVertex:mode = 2, form_time records the late formation time from
 * Particle::tProd(), and prod_time records the earlier breakup time from
 * Particle::tProduction(), including origin-only production assignments.
 * The appended prod_time_status column is 1 for an available production
 * reference (including zero), 0 for unavailable information (prod_time = 0),
 * and 2 for not applicable (spectator nucleons, with both times zero).
 * Non-finite production references are errors rather than status values.
 *
 * Units: input particle vertices in mm and coordinate times in mm/c;
 * output positions in fm and t, prod_time and form_time in fm/c. Momentum,
 * energy and mass remain in GeV. Coordinate times refer to the input particle
 * frame and are distinct from proper-time intervals. The complete saved vertex
 * and the separate production time are converted to fm before freestreaming.
 * The writer adds no spatial smearing to hadron positions.
 *
 * @param event_number Global event ID.
 * @param atomic_number Atomic number \f$ Z \f$ of the target nucleus.
 * @param mass_number Mass number \f$ A \f$ of the target nucleus.
 * @param struck_nucleon Nucleon assigned by the hadronizer remnant model.
 * @param kin Reconstructed DIS kinematics for the event.
 * @param particles PYTHIA particles returned by the hadronizer; vertex
 *        units are mm, with coordinate times in mm/c.
 * @param event_out Stream that receives the OSCAR2013 particle list.
 * @param meta_out Stream that receives the JSON metadata record.
 * @param history Optional observer of the same output positions and spectators;
 *        consumes no additional random draws.
 * @param target Optional shared target for deterministic untouched spectators.
 * @param spectator_rng Persistent seeded engine, required when legacy spectators
 *        are emitted. Explicit targets and events without spectators need none.
 * @param reproducibility Optional JSON object containing seed/input provenance.
 */
void write_event_output(
    int32_t event_number,
    int atomic_number,
    int mass_number,
    StruckNucleon struck_nucleon,
    const DISKinematics& kin,
    const std::vector<Pythia8::Particle>& particles,
    std::ostream& event_out,
    std::ostream& meta_out,
    SpacetimeHistory* history = nullptr,
    const EHIJING::NuclearTarget* target = nullptr,
    std::mt19937* spectator_rng = nullptr,
    const std::string& reproducibility = "");

// Actual untouched identities must agree with final B/Q, without deleting or
// inventing spectators. A mismatch rejects the whole explicit-target event.
bool target_balance_matches(const EHIJING::NuclearTarget& target,
                            const SpectatorBalance& balance);
