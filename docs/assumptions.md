# Physical assumptions and open questions {#physics_assumptions}

This page records the assumptions, known limitations and unresolved consistency
questions in the eHIJING checkout used by ELECTRA. It complements the
[physics manual and execution diagram](mainpage.md) and the
[unit and frame conventions](units.md). The theoretical background is the
[eHIJING paper](https://arxiv.org/html/2304.10779v2); agreement with that paper's
model is not, by itself, validation of subsequent ELECTRA extensions.

**Review basis: 14 October 2026.** Statements below are based on source,
settings, existing documentation and synthetic tests. The C4 resolution also
includes executed fresh-process replay checks with small temporary eHIJING runs
and synthetic SMASH transport. These checks establish reproducibility in their
stated scope; they are not production-level physics benchmarks. Proposed checks
in the open items remain work to be done. This page chooses no new model parameters.

Paths beginning with `src/`, `include/`, `tests/` or `external/` are relative
to the eHIJING component. Paths explicitly described as ELECTRA paths are
relative to the parent ELECTRA repository.

## How to interpret the status

| Status | Meaning |
| --- | --- |
| Model assumption | An intentional physical approximation or prescription. Its usefulness has a domain of applicability; it is not automatically a software defect. |
| Confirmed limitation | A behavior or missing guarantee directly established by the current source. Its size or phenomenological impact may still be unmeasured. |
| Confirmed mismatch | Two current source interfaces or normalizations disagree under the stated conditions. This does not establish which installed binary or configuration was used in an earlier run. |
| Resolved configuration mismatch | The checked-in configuration now matches the source interface. Runtime and physics validation remain separate checks. |
| Resolved for fixed chunks | Seeded streams and complete-record replay are established for fixed chunk starts with identical inputs, build and environment; independent event replay after rechunking is outside the guarantee. |
| Open consistency question | There is a specific reason to check consistency, but this review has not established the outcome with an independent calculation or runtime test. |
| Validation gap | Existing tests do not establish the stated physical claim. Absence of that validation is not proof that the result is wrong. |

An item should be closed only with evidence that addresses its actual scope.
For example, conserving momentum in an isolated elastic collision does not
establish conservation in a radiative event, and converting mm to fm does not
establish frame consistency.

## Physical assumptions currently made

### Hard process and the nuclear initial state

The hard event is generated for a per-nucleon target using PYTHIA. In the
supplied HERMES setup it is neutral-current DIS, with nuclear PDFs or isospin
corrections supplied through the hard-process PDF interface. The driver selects
EPPS16 NLO for `Z > 2` and isospin only otherwise; the hard matrix element does
not thereby become an NLO calculation.

The primary nuclear identity is a geometrical assignment. Its species is
selected with probability `Z/A`, then an identity is selected within that
species. Its center supplies the hard vertex, without resampling the hard
PDF, assigning a Fermi momentum or accounting for nuclear binding energy.
Consequently, fitted nPDF effects and explicit nucleon geometry are separate
ingredients rather than one jointly sampled nuclear wave function.

The localized hard-vertex picture and predominantly external fragmentation
are motivated by the approximate hierarchies

\f[
 \frac{1}{Q}\ll L,\qquad
 \frac{1}{2m_Nx_B}\ll L,\qquad
 \frac{z_h\nu}{\Lambda_{\rm QCD}^2}\gg L.
\f]

These conditions are not automatically imposed by the event-selection code.
At small Bjorken x there is no switch to a coherent photon-dipole interaction
with the full nucleus. The current hadronic charge budget supports
neutral-current events; charged-current bookkeeping needs a separate treatment.

**Implementation:** `src/main.cpp`, `src/dis_kinematics.cpp`,
`external/pythia8/src/PartonLevel.cc` and
`external/pythia8/include/eHIJING/NuclearTarget.h`.

### Frozen geometry and eikonal trajectories

Explicit-target mode samples exactly A nucleons, with fixed proton/neutron
identities. Heavy-target centers follow the configured spherical Woods--Saxon
prescription; deuterium uses the Hulthén separation prescription. Gaussian
profiles of width `--nucleon-width-fm` supply the density for rescattering.
The default width of 0.5 fm is a benchmark, not a fitted uncertainty estimate.

The primary identity is excluded from medium scattering. Other centers are
independent: there is no short-range exclusion, recentering, nuclear
deformation, Fermi-motion or binding model. The density remains frozen when
the remnant momenta change. Collision locations are sampled on straight
forward rays, and kicks do not cause those already sampled paths to be
recomputed. In particular, exact elastic momentum recoil does not introduce
a dynamically moving nuclear density.

Legacy `--nuclear-target off` uses a different continuous hard-sphere
prescription, including a full-chord thickness, its subsequent rescaling and
an initial collision-time exclusion interval. It is an alternative
prescription, rather than an exact limit of the Gaussian identity model.

**Implementation:** `NuclearTarget::sample()`, `path_weights()` and
`sample_collision()`, plus `EHIJING::NuclearGeometry` and the two collision
samplers.

### Microscopic collisions and medium radiation

The medium rate uses a saturation-motivated gluon TMD with parameters
K, n and lambda. It is a phenomenological input without dynamical target
small-x evolution. A finite transfer cutoff resolves a finite number of
collisions. Independent eikonal sampling is subsequently combined with
collision-conditioned radiation weights.

HT and GHT are alternative medium kernels in the soft-gluon approximation.
They do not constitute two independent corrections to be added together.
Finite-z numerators and factors of `1-z` in the implementation do not establish
a full finite-z derivation. The low-scale adapter modifies light quarks,
antiquarks and gluons; it is not a dedicated heavy-quark transport or
heavy-quark medium-radiation model.

The saturation scale also supplies a high/low emission separation
prescription. The explicit-target shower floors that separation at `mu2`,
while the saturation solver itself can return smaller roots. Below the
separation scale, the adapter uses a frozen coupling evaluated at that scale
and an energy resolution of 0.2 GeV. With the implemented freeze prescription,
the limiting medium coupling is about 1.40. Thus the cutoff region should not
be interpreted as a controlled weak-coupling expansion simply because its
kernel has a perturbative origin.

**Implementation:** `external/pythia8/src/eHIJING.cc`,
`external/pythia8/src/TimeShower.cc` and `src/modified_ff.cpp`.

### Remnants, hadronization and downstream transport

Explicit-target recoil treats each touched nonprimary nucleon as a shared
on-shell remnant. After its hits, it is split once into an isotropic
constituent quark--diquark pair with the fixed species' valence choices.
Color exchanges reuse those remnant strings. Untouched identities become
spectators. There is no model of the residual nucleus' binding, evaporation,
excitation cascade or collective motion.

Hadronization uses a separate PYTHIA instance and Lund fragmentation with the
selected tune and settings. Production references and late formation vertices
are assignments of that string model. They are not direct measurements of
when a color-neutral object first interacts with the nucleus. SMASH then adds
its configured hadronic interaction and formation-profile model; that is a
further physical prescription, not a consequence of unit conversion.

**Implementation:** `src/modified_ff.cpp`, `include/nuclear_recoil.hpp`,
`src/hadronizer.cpp`, the bundled string fragmentation sources, and
`src/event_output.cpp`. The transport interface is in ELECTRA's
`src/smash/src/listmodus.cc`.

## Review register

The identifiers below are intended to remain stable as evidence is added.
The order separates source-confirmed limitations from questions still requiring
an independent check; it is not a measured ranking of their numerical effects.

| ID | Status | Subject |
| --- | --- | --- |
| C1 | Confirmed limitation | Low-scale radiation does not guarantee four-momentum closure. |
| C2 | Confirmed limitation | High-scale shower delays are lost during hadronizer import. |
| C3 | Confirmed mismatch | The two legacy qhat helpers have different dimensions. |
| C4 | Resolved for fixed chunks | A versioned master seed controls all active event streams; rechunking is outside the guarantee. |
| C5 | Resolved configuration mismatch | The checked-in SMASH template now declares the production and formation columns required by the current reader. |
| C6 | Resolved output ambiguity; explicit fallback | Production-time status distinguishes missing information, available times and spectators; SMASH warns when using the zero fallback. |
| C7 | Confirmed limitation | Cached table files do not establish their parameter provenance. |
| Q1 | Open consistency question | How much does rejection condition the accepted physics sample? |
| Q2 | Open consistency question | Does the transfer proposal/veto sequence reproduce the intended rate and cutoff? |
| Q3 | Open consistency question | Is the shared separation scale consistent across radiator paths and shower stages? |
| Q4 | Open consistency question | Are the GHT approximations and veto bounds reliable across their full domain? |
| Q5 | Open consistency question | Are combined formation assignments and frame interfaces physically consistent? |
| Q6 | Open consistency question | Is primary species assignment consistent with the hard PDF's flavor mixture? |
| Q7 | Open consistency question | Is the PYTHIA recoil baseline suitable for the medium modifications used here? |
| V1 | Validation gap | The modified generator and its SMASH handoff need physics benchmarks beyond synthetic tests. |

## Confirmed limitations and mismatches

### C1. Radiation and global energy--momentum closure

**Evidence.** In `Modified_FF::sample_ff_partons()`, a sampled gluon's spatial
momentum is subtracted from the parent, whose energy is then recomputed
independently on shell:

```cpp
particle.p(particle.p() - kmu);
particle.e(std::sqrt(particle.pAbs2() + particle.m2()));
kmu.e(std::sqrt(kmu.pAbs2()));
```

Consequently, even when the parent and gluon are individually on shell,

\f[
 \Delta E_{\rm rad}
 =\sqrt{|\mathbf p-\mathbf p_g|^2+m^2}+E_g-E
\f]

is generally nonzero. `recoil_two_body()` conserves the four-momentum and
masses of each subsequent accepted elastic collision, but it does not repair
this radiative residual. The final baryon-number and charge checks do not
test energy--momentum closure.

**Physical consequence.** Total energy, reconstructed hadronic masses and
fragmentation observables can inherit an error whose size has not been
established in this review. It is therefore inappropriate to claim exact
global event conservation from the elastic solver's tests alone.

**Evidence needed to resolve it.** Record the four-momentum residual at each
stage with consistent incoming/outgoing particle sets and a common frame.
Separate radiation, elastic recoil, remnant splitting, fragmentation and
output contributions. Quantify their distributions before selecting any
new radiative recoil prescription.

### C2. Missing high-scale shower spacetime history

**Evidence.** `Hadronizer::hadronize()` imports final partons with status 23,
momenta, colors and stored displacements, but no source parent indices.
Its ancestry traversal therefore cannot reconstruct high-scale shower delays;
the corresponding accumulated displacement remains zero. Existing medium
displacements and subsequently generated Lund vertices survive.

**Physical consequence.** A formation position in the output is not the sum
of a reconstructed hard-shower worldline, medium transport and string
fragmentation. Treating the present output as that complete spacetime history
would overstate the information available. HepMC observation snapshots do not
restore the missing operative ancestry.

**Evidence needed to resolve it.** Specify the intended shower-lifetime
prescription and its reference frame, preserve the necessary ancestry, and
verify that hard-origin translations and medium displacements are each added
once. This requires a physical modeling decision as well as bookkeeping.
Unit conversion alone cannot close this item.

### C3. Transport-coefficient helper normalization

**Evidence.** `MultipleCollision::qhatA()` returns
`rho0 * Qs2 / TA`, whereas `qhatF()` returns
`(CF/CA) * Qs2 / TA`. With density in GeV cubed and thickness and saturation
scale in GeV squared, the former has transport-coefficient dimensions and the
latter is dimensionless.

**Physical consequence.** A result obtained by calling the latter helper
cannot be labeled GeV squared per fm without supplying the missing
normalization. This observation does not establish that the helper is used
in an existing physics analysis; its downstream consumers must be traced.

**Evidence needed to resolve it.** Audit its callers and compare reported
values with the rate definition and color-factor relation,

\f[
 \hat q_R=\int d^2\boldsymbol\kappa\,\kappa^2
             \frac{d\Gamma_R}{d^2\boldsymbol\kappa},\qquad
 \frac{\hat q_F}{\hat q_A}=\frac{C_F}{C_A},
\f]

using identical cutoffs and paths. Keep collisional and radiative transport
coefficients distinct, and state the conversion from natural units to
GeV squared per fm.

### C4. Seed propagation and reproducibility

**Resolution.** `electra-seeds-v1` preserves a validated unsigned 64-bit
master seed at the CLI boundary. The executable derives each stream from that
master, the chunk's first event ID and a frozen stream identifier. PYTHIA
receives integer seeds in `[1, 900000000]`; all medium, hadronizer-helper,
colour-choice and legacy-spectator engines receive explicit seeds. The explicit
nuclear target continues to use the hard PYTHIA stream. No event stream obtains
independent entropy. Omitting the standalone seed obtains entropy once for the
master, which is printed and recorded.

**Policy across attempts, events and chunks.** Streams are initialized before
sampling and advance continuously through all failed attempts and accepted
events. They are never reset on rejection. The successful attempt number is
monotonic for each requested event, including across the driver's 15-attempt
batches. Event count and output shard width are excluded from derivation, so
extending a fresh run with the same chunk start preserves its prefix. A
different chunk start changes its stream context. This resolves repeated-run
reproducibility for fixed chunk boundaries; it does not make events independently
replayable after rechunking or changing physics parameters.

**Evidence.** `tests/seed_policy_test.cpp` checks fixed derivation vectors,
wide-seed parsing, engine replay and PYTHIA's seed boundaries. The parent
repository's `tests/reproducibility_smoke.py` compares full OSCAR, metadata,
fragmentation and HepMC3 records across fresh processes for explicit/legacy
targets and HT/GHT modes. It also checks diagnostic collection, direct/wrapper
agreement, changed seeds and prefix preservation. The parent repository's
`tests/smash_seed_smoke.py` also compares complete synthetic transport records and
both replica seeds across
fresh processes and task modes, and verifies that changing the master changes
the stochastic records. Qs/GHT tables now serialize with round-trip precision;
generation uses locking and atomic publication so
parallel chunks do not load partial tables. Existing tables are preserved.

**Scope and provenance.** Identical settings, table bytes, executable and
runtime environment are required. Metadata includes the policy, master and
effective stream seeds, successful attempt and input SHA-256 digests. Completed
wrapper tasks validate saved requests before skipping; historical completion
markers without provenance require a new run directory. Cross-architecture
bitwise equality, disjoint random streams and matched parameter-variation
histories are not claimed. See the [reproducibility policy](../README.md#reproducibility-and-seeds)
for the stream inventory and exact derivation.

### C5. SMASH formation-field configuration mismatch

**Evidence.** The current modified reader in ELECTRA's
`src/smash/src/listmodus.cc` requires both `prod_time` and `form_time` in
`Modi:List:Optional_Quantities`, requires values on each row, and passes them
directly to `set_slow_formation_times()`. The eHIJING writer exports both
columns. ELECTRA's checked-in `input/smash/config_files/config.yaml` now declares

```yaml
Optional_Quantities: ["ID", "charge", "prod_time", "form_time", "xsecfac", "prod_time_status"]
```

The Python `smash_physical_event_task()` overrides the input filename,
random seed, event count and profile path, but does not replace that optional
field list. `slurm-scripts/smash_submit.sh` defaults to this template.

**Resolution.** The obsolete `begin_form_time` entry has been replaced with
`prod_time` and `form_time`, in the same order as the eHIJING writer's OSCAR
header and particle rows. This resolves the checked-in source/configuration
mismatch without changing the formation model or its parameters. The profile
`Formation_Time` parameter remains configured, but does not replace or rescale
the coordinate times required by list input.

**Remaining validation.** No SMASH run was attempted for this configuration
fix. Check the exact installed reader binary, universal configuration and
effective command-line overrides together, then perform a small interface
round trip before a physics production. Custom configurations or older
containers may use another interface.

### C6. Missing versus zero production times

**Original limitation.** `write_final_hadrons()` in `src/event_output.cpp`
wrote zero when `hasProductionTime()` was false without retaining availability
in the OSCAR row. An assigned zero and unavailable production information were
therefore indistinguishable to the reader, although diagnostics retained the
production-reference flag.

**Resolution.** Every OSCAR particle row now appends the dimensionless integer
`prod_time_status` after `xsecfac`:

| Value | Meaning | Exported production time |
| --- | --- | --- |
| `0` | Production reference unavailable | Zero placeholder |
| `1` | Production reference available, including origin-only assignments | Stored coordinate time, which may be zero |
| `2` | Not applicable: pre-existing spectator nucleon | Zero, with formation time also zero |

The status comes from the final exported particle's `hasProductionTime()`
flag, with spectators handled separately in both legacy and explicit-target
output. Non-finite available production references cause an error rather than
an exceptional status. The optional history retains the same value as
`oscar_prod_time_status`.

The checked-in SMASH template declares the new column in `Optional_Quantities`.
`src/smash/src/listmodus.cc` validates status values and forces the numerical
profile start time to zero for status `0`, even if a nonzero placeholder was
supplied. It prints a warning to stderr identifying the PDG species and the
zero fallback. Status `1` and `2` use the supplied numerical times without
that warning. Warnings occur during startup validation and particle import.
Legacy files configured without the status column retain their numerical
times, without inferring availability or emitting this warning. Rebuild the
writer and reader binaries/container and update custom configurations to use
the new schema; configuration order must match the particle rows.

**Validation.** Synthetic writer/history checks cover missing references,
available zero and nonzero times, origin-only assignments, both spectator
paths and non-finite references. Focused SMASH reader tests cover the new
particle-row schema, field order, stderr warnings, zero fallback, malformed or
missing statuses, and legacy numerical-time input. No simulation events were
generated for this interface change.

**Physical consequence and remaining scope.** The missing-versus-zero output
ambiguity is resolved for the new schema. The zero fallback is an explicit
transport policy, not evidence of production at the hard vertex. This change
does not reconstruct missing references or establish the physical production
event of every hadron or decay product. Missing-reference rates by fragmentation
branch, copy/decay history and final species remain a separate physics audit.

### C7. Table provenance, saturation roots and scale floors

**Evidence.** The table loader reads existing `Qs.dat` and `GHT.dat` files
using their expected sizes. `Qs.dat` does not encode or validate the K, n and
lambda values with which it was generated. The direct solver can produce
roots below `mu2`, while the explicit-target high-scale separation is floored
at `mu2`; out-of-domain explicit-target queries are solved directly.

**Physical consequence.** Reusing an unrelated table directory can make the
cached saturation input inconsistent with the requested TMD parameters. A
valid table size does not establish physical provenance. The raw saturation
root, the elastic screening scale and the stored shower separation also need
to be distinguished when interpreting dilute or surface trajectories.

**Evidence needed to resolve it.** Retain parameter values, source revision,
grid limits and a content hash with each table set. Compare interpolation
against direct roots near grid boundaries and sparse-density regions, and
check continuity where the matching floor becomes active. This review has
not established that a particular existing table set is stale.

## Open consistency questions

### Q1. Selection bias from rejecting entire attempts

**Evidence and uncertainty.** An impossible forward elastic solution rejects
the whole attempt. Color/target-state failures, hadronization failures and
final target-ledger mismatches can also trigger regeneration. The effect on
the accepted kinematic and geometric distribution has not been measured here.
If the proposed distribution is g and its acceptance probability is a, then

\f[
 p_{\rm accepted}(\xi)
 =\frac{g(\xi)\,a(\xi)}{\int d\xi'\,g(\xi')\,a(\xi')}.
\f]

Thus a low average rejection rate alone does not establish an unbiased result
in a particular phase-space bin.

**Check needed.** Resolve rejection reasons versus DIS kinematics, parton
energy, transfer, target species, thickness, hit multiplicity and geometry.
Compare proposed and accepted distributions. Determine whether the resulting
conditioning is an intended physical phase-space restriction or exposes a
mismatch between the proposal kernel and the exact recoil prescription.

**Sources:** `src/main.cpp`, `src/modified_ff.cpp` and
`include/nuclear_recoil.hpp`.

### Q2. Transfer proposal, veto and cutoff consistency

**Evidence and uncertainty.** In both collision samplers, `xg` is computed
from the previous value of `q2`, then `q2` is updated and the veto uses the
stored `xg`. The loop condition tests the previous transfer against `q2min`;
there is no new explicit lower-cutoff check between the update and append.
For example, the explicit-target sampler contains

```cpp
const double xg = q2 / q2max;
q2 *= std::pow(1. + (lambdaG_ - 1.) * lnr * q2 / tildeTA
    * std::pow(xg, -lambdaG_), 1. / (lambdaG_ - 1.));
// Keep the original eHIJING veto kernel and momentum-transfer proposal.
if (flat_gen(gen) < std::pow(1. - xg, powerG_) * q2 / (q2 + qs2)) {
```

These source observations justify an independent sampling audit. This page
does not claim that a distributional discrepancy has already been quantified,
or silently replace the algorithm with another interpretation.

**Check needed.** Derive the precise proposal and veto measure using old and
new transfer variables explicitly. Compare sampled spectra and multiplicities
with direct quadrature of the intended rate on a fixed medium, and test
whether accepted transfers remain above the advertised resolution. Test the
zero-medium limit and the endpoints separately from recoil acceptance.

**Sources:** `MultipleCollision::sample_all_qt2()` and
`sample_all_qt2_target()` in `external/pythia8/src/eHIJING.cc`.

### Q3. Separation scale and the ordering of the two shower stages

**Evidence and uncertainty.** `TimeShower` computes a ray-dependent scale
for the current radiator and updates the event's single `SeparationScale()`.
`Modified_FF` subsequently reads that stored value once and uses it across
its final-parton loop. Its radiation is sampled before the collected elastic
updates are sorted and applied. A newly emitted gluon's initial displacement
comes from a formation ordering variable, and its extra collisions begin
there. This does not establish a unique, coupled time evolution for the event.

**Check needed.** Record which radiator last set the stored scale and compare
it with the final partons' path-dependent scales. Check for stage-overlap or
coverage gaps using the actual `kt2`, `lt2` and PYTHIA ordering variables.
Establish which corrections are valid when radiation and recoil change the
parent momentum used by subsequent operations. A different matching or
time-evolution prescription would require an explicit model choice.

**Sources:** `TimeShower::pT2nextQCD()` and
`Modified_FF::sample_ff_partons()`.

### Q4. GHT angular approximations and probability bounds

**Evidence and uncertainty.** The high-scale GHT implementation switches
between a Bessel-function approximation, angular-table interpolation and
an asymptotic expression. Its returned weight is clipped to the interval
from zero to one. The low-scale sampler also uses ratios of interference
weights and supporting-collision sums. The existence of those formulas and
clamps does not establish that all overestimates remain valid or that limiting
cases are handled accurately. An angular interference contribution need not
itself be positive.

**Check needed.** Compare each approximation with direct angular integration,
including branch boundaries, tiny formation phases, nearly equal radiation
and transfer scales, and endpoint energy fractions. Inspect un-clipped
weights and count invalid, nonfinite or out-of-range trials. Verify that
the proposal dominates the intended probability and that any clipping has
negligible, understood impact before treating it as harmless numerics.

**Sources:** `compute_GHT_Angular_Table()`, `sample_z_stochastic()` and the
low-scale HT/GHT acceptance branches in `src/modified_ff.cpp`.

### Q5. Combined spacetime assignments and frame consistency

**Evidence and uncertainty.** Explicit-target interfaces boost full momenta
and four-positions, and synthetic tests check selected boosted recoil
displacements. However, shower formation variables are ordering variables,
the hard-shower delays in C2 are absent, and Lund production references can
be inherited origins or breakup vertices depending on the string branch.
A collapse's inherited average origin is not a reconstructed color-exchange
event. These ingredients do not yet establish a complete causal spacetime
interpretation for every exported particle.

The writer converts the saved vertices to fm and extrapolates to t = 0.
The modified SMASH reader uses the imported production and formation values
as coordinate times, without another boost factor. Correct conversions are
necessary but do not prove that every contributing assignment shares the
intended physical frame and origin.

**Check needed.** For each fragmentation branch, compare the full production
references, saved formation vertex, inherited origin and the output surface
in a common frame. Inspect nonfinite times and the ordering of profile
endpoints. Use controlled boosts to distinguish covariance of interface
four-vectors from the separate question of equivalence of the full sampling
algorithm. Quantify how these assignments affect the downstream formation
profile and possible overlap of partonic and hadronic interaction periods.

**Sources:** `src/hadronizer.cpp`, `src/event_output.cpp`, bundled
`Event.h`, `StringFragmentation.cc`, `MiniStringFragmentation.cc`,
`docs/units.md`, and ELECTRA's modified SMASH list reader.

### Q6. Primary identity and the PDF flavor mixture

**Evidence and uncertainty.** The geometry chooses the primary species with
probability Z/A independently of the generated hard flavor. The hadronizer
reuses that species for primary remnant construction. This geometrical choice
is not explicitly conditioned on the proton and neutron contributions to
the hard PDF for the selected flavor and kinematics.

**Check needed.** Compare the joint distribution of primary species and hard
flavor with the intended per-nucleon nPDF mixture, including valence and sea
channels. Check the final B/Q acceptance separately by primary species and
flavor. Agreement of the unconditional proton fraction with Z/A, or of final
net B/Q alone, does not establish that conditional flavor correlations are
correct. The physical impact remains to be determined.

**Sources:** `NuclearTarget::sample()`, hard-PDF setup in `src/main.cpp`,
and primary-remnant handling in `src/hadronizer.cpp`.

### Q7. DIS recoil scheme and the vacuum reference

**Evidence and uncertainty.** The current HERMES hard-process settings enable
`SpaceShower:dipoleRecoil` and disable ISR/FSR interleaving. The modified
`TimeShower` retains FSR when eHIJING is active. The physics manual already
distinguishes these settings from the paper's global-recoil setup. A model
comparison using the paper's parameters does not establish that this modified
vacuum/medium combination is equivalent.

**Check needed.** Establish the active recoil and shower settings for the
actual linked PYTHIA version. Validate the vacuum DIS baseline, reconstructed
lepton kinematics and fragmentation before comparing medium effects.
Quantify recoil-scheme sensitivity under an explicitly selected comparison
plan; do not assume that changing one setting restores the published setup.

**Sources:** ELECTRA's
`input/ehijing/experiments/hermes/hard_process.setting`, the modified
`TimeShower.cc`, and `src/dis_kinematics.cpp`.

## What existing tests establish

The table describes coverage in the test sources. The C4 replay checks were also
executed for the seed-policy resolution, as recorded in that item; their temporary
fixtures do not replace physics validation with a production sample.

| Test source | Coverage present | What that coverage does not establish |
| --- | --- | --- |
| `tests/nuclear_target_test.cpp` | Composition and identity rules, selected Gaussian integrals and collision-position statistics, interpolation endpoint, elastic masses and four-momentum in constructed cases, impossible-transfer rejection. | The full nonzero-medium transfer distribution, a radiative shower, or measured nuclear-observable agreement. |
| `tests/shared_remnant_test.cpp` | Repeated hits share one remnant pair, chronological elastic updates, a constructed boosted displacement, conservation in the nonradiative fixture, and rejection of an impossible transfer. | Radiation energy closure: its soft-jet fixture deliberately avoids radiation and uses an artificial saturation table. |
| `tests/spectator_balance_test.cpp` | Signed final B/Q counting, resonances and antiparticles, explicit untouched identities, output counts, and labeled legacy fallback behavior. | Exact total energy, full flavor correlations, nuclear binding or an unbiased accepted event sample. |
| `tests/spacetime_history_test.cpp` | HepMC round trips, record identity, available versus missing vertices, coordinate conversions and diagnostic output invariance for synthetic fixtures. | Reconstruction of absent shower lifetimes, causal ordering of all assigned vertices, or the physics of a SMASH formation profile. |
| `tests/seed_policy_test.cpp` | Frozen derivation vectors, unsigned 64-bit parsing, engine replay, high seed bits, PYTHIA boundaries and SHA-256 reference vectors. | Disjoint streams, cross-platform event equality or physical validity of the generated sample. |
| ELECTRA `tests/test_seed_policy.py` | Master validation, unambiguous structured hashing, unchanged forwarding, recorded effective seeds, completion fingerprints and captured SLURM submission commands. | A real SLURM/container run or physics validation; subprocess execution is mocked in the forwarding checks. |
| ELECTRA `tests/reproducibility_smoke.py` | Complete-record replay in fresh processes for HERMES deuterium fixtures with HT/GHT and explicit/legacy targets; table generation/reload, observers on/off, path aliases, prefix preservation and direct/wrapper agreement. | Rechunking independence, cross-architecture equality, exhaustive coverage of target/parameter settings or phenomenological agreement. |
| ELECTRA `tests/smash_seed_smoke.py` | Complete synthetic pion transport records and two replica seeds agree across fresh processes/task modes; a changed master changes stochastic output. | Production transport benchmarks, arbitrary profiles/inputs or the physical interpretation of formation times. |

### V1. Production-level physics validation remains open

The component README records that the explicit-target extensions have not
yet been validated with a physics production sample. The tests above provide
useful component checks, but none closes that gap. In particular, a K value
used with the published geometry and shower setup is not automatically a
validated parameter choice for the modified geometry, recoil and downstream
transport combination.

A focused validation sequence would be:

1. Establish the exact source/library/settings versions, table provenance,
   effective SMASH fields and random-stream policy. Resolve interface
   mismatches before interpreting physical output.
2. Check elementary rates and samplers independently: transfer distributions,
   saturation roots, HT/GHT weights, probability bounds and dimensional
   normalizations.
3. Extend controlled fixtures to combined radiation/recoil, complete stage
   four-momentum accounting, and production/formation references through
   boosts, copies, fragmentation branches and the output/reader boundary.
4. Under an agreed physics comparison plan, compare proposed and accepted
   distributions, the vacuum reference, geometry changes and the two
   medium kernels with matched cuts and recorded rejection reasons.
5. Compare physical observables with appropriate reference data and report
   both parameter and model dependence. Separate eHIJING-only results from
   changes introduced by SMASH and its formation profile.

These are proposed investigations, not results obtained for this page.
Changing geometry, recoil, matching, fragmentation or formation prescriptions
requires choosing the physical model and documenting its alternatives before
implementation.

## Maintaining this review

For each resolved item, retain the affected revision and configuration,
independent calculation or test, numerical result, and the scope of the
conclusion. A useful record distinguishes an implementation repair from a
new model choice and from a phenomenological fit. Keep the stable identifier
and state which remaining observable or frame range is still unvalidated.
Do not close a physics question solely because documentation builds or a
serialization test pass.
