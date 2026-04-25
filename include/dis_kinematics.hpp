#pragma once

#include "Pythia8/Pythia.h"

/**
 * @brief Container for reconstructed deep-inelastic scattering kinematics.
 *
 * Stores the incoming and outgoing four-momenta needed to reconstruct the
 * exchanged boson momentum and the standard Lorentz-invariant DIS variables.
 *
 * The notation follows the usual DIS convention:
 *
 * - \f$ P \f$: incoming target/proton four-momentum
 * - \f$ k \f$: incoming lepton four-momentum
 * - \f$ k' \f$: outgoing scattered lepton four-momentum
 * - \f$ q = k - k' \f$: exchanged virtual boson four-momentum
 */
struct DISKinematics {
    /// Incoming proton/target four-momentum, \f$ P \f$.
    Pythia8::Vec4 pProton;

    /// Incoming lepton four-momentum, \f$ k \f$.
    Pythia8::Vec4 peIn;

    /// Outgoing scattered lepton four-momentum, \f$ k' \f$.
    Pythia8::Vec4 peOut;

    /// Exchanged virtual boson four-momentum, \f$ q = k - k' \f$.
    Pythia8::Vec4 pGamma;

    /// Energy transfer in the target rest frame, \f$ \nu = P \cdot q / \sqrt{P^2} \f$.
    double nu = 0.0;

    /// Photon virtuality, \f$ Q^2 = -q^2 \f$, in \f$ \mathrm{GeV}^2 \f$.
    double Q2 = 0.0;

    /// Positive square root of the photon virtuality, \f$ Q = \sqrt{Q^2} \f$, in GeV.
    double Q = 0.0;

    /// Invariant mass squared of the photon-target system, \f$ W^2 = (P + q)^2 \f$.
    double W2 = 0.0;

    /// Invariant mass of the photon-target system, \f$ W = \sqrt{W^2} \f$.
    double W = 0.0;

    /// Bjorken scaling variable, \f$ x_B = Q^2 / (2 P \cdot q) \f$.
    double bjorken_x = 0.0;

    /// DIS inelasticity, \f$ y = (P \cdot q)/(P \cdot k) \f$.
    double y = 0.0;
};

/**
 * @brief Compute DIS kinematics from the current PYTHIA event record.
 *
 * Extracts the incoming proton, incoming lepton, and outgoing lepton from the
 * PYTHIA event record, constructs the exchanged boson momentum,
 * \f$ q = k - k' \f$, and computes the standard DIS variables
 * \f$ \nu \f$, \f$ Q^2 \f$, \f$ Q \f$, \f$ W^2 \f$, \f$ W \f$,
 * \f$ x_B \f$, and \f$ y \f$.
 *
 * @param pythia PYTHIA generator instance containing the current event.
 *
 * @return A DISKinematics object containing the relevant four-momenta and
 *         derived DIS variables.
 *
 * @note This function assumes the following PYTHIA event-record convention:
 *       - `event[1]`: incoming proton/target
 *       - `event[4]`: incoming lepton
 *       - `event[6]`: outgoing scattered lepton
 *
 * @warning If the event record does not follow this indexing convention, the
 *          reconstructed DIS kinematics will be incorrect.
 */
DISKinematics compute_dis_kinematics(const Pythia8::Pythia& pythia);

/**
 * @brief Check whether reconstructed DIS kinematics pass the event selection.
 *
 * Applies the DIS phase-space cuts used by the analysis/generator filter.
 *
 * @param kin Reconstructed DIS kinematics for a single event.
 *
 * @return `true` if the event passes all selection cuts; otherwise `false`.
 */
bool is_valid_dis_event(const DISKinematics& kin);

/**
 * @brief Apply the DIS event trigger to the current PYTHIA event.
 *
 * Convenience wrapper that computes the DIS kinematics from the current PYTHIA
 * event and immediately checks whether the event passes the DIS selection cuts.
 *
 * @param pythia PYTHIA generator instance containing the current event.
 *
 * @return `true` if the current event passes the DIS selection; otherwise
 *         `false`.
 *
 * @see compute_dis_kinematics
 * @see is_valid_dis_event
 */
bool trigger(const Pythia8::Pythia& pythia);