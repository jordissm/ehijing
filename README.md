## Documentation

The eHIJING headers and sources are documented with Doxygen.

[Learning the physics of eHIJING](docs/physics-tutorial.md) is a graduate-level
tutorial that starts with DIS and QCD showers, derives the medium ingredients,
and connects them to this implementation through worked examples, code excerpts
and exercises.

[Physics included in eHIJING](docs/mainpage.md#physics-included-in-ehijing)
details DIS kinematics, nuclear PDFs and geometry, the gluon TMD and saturation
scale, multiple scattering, HT/GHT radiation and LPM interference, shower
matching, recoil, and Lund hadronization, with LaTeX equations and an explicit
account of the implementation's approximations.

[Physical assumptions and open questions](docs/assumptions.md) separates
model choices, source-confirmed limitations, unresolved consistency questions
and the validation still needed for the current implementation.

[Unit conventions](docs/units.md) describe spacetime units, conversion factors,
interface requirements, and remaining review items.

## Shared nuclear target and recoil

The ELECTRA executable now defaults to `--nuclear-target on`. Rebuild both
`external/pythia8` and the CMake application, and install the updated PYTHIA XML
settings together with the library. A link-time ABI guard prevents linking the
application to an old library. The underlying PYTHIA setting
`eHIJING:ExplicitNucleus` defaults off for other library clients. Use
`--nuclear-target off` to retain ELECTRA's previous continuous-medium recoil
and spectator-count prescription.

Each generation attempt constructs exactly A identified nucleons, with exactly
Z protons. The primary species is chosen with Z/A, then an identity is chosen
uniformly within that species. Its centre is the hard vertex, and its species
is reused by the hadronizer. The primary identity is excluded from subsequent
medium scattering. Hard-process momenta and nuclear PDFs retain their existing
sampling; this is a geometrical assignment, without a Fermi-motion model.

The geometry follows [TRENTo's nucleus prescriptions](https://github.com/Duke-QCD/trento/blob/master/src/nucleus.cxx)
and its [Gaussian nucleon profile](https://github.com/Duke-QCD/trento/blob/master/doc/usage.rst):
A=1 is a single nucleon; A=2, Z=1 uses the Hulthén proton-neutron separation
with a=0.228 and b=1.18 fm^-1, placing nucleons at opposite half-separations.
Heavy targets use spherical Woods–Saxon centres with the volume Jacobian and
radial cutoff R+10a. Centres are independent, without minimum-distance
correlations or recentering. Density is the sum of normalized spherical
Gaussians with `--nucleon-width-fm 0.5` (a benchmark, not a fitted value).
Geometry remains frozen; recoil momenta do not move the density centres.
Collision sampling retains straight eikonal rays; kicks do not resample paths. No
hotspots, strength fluctuations, nuclear binding, deformations, or TRENTo3D
heavy-ion energy-deposition prescription are added.

Default Woods–Saxon parameters in fm are:

| A | R | a |
| --- | --- | --- |
| 63 | 4.20 | 0.596 |
| 129 | 5.36 | 0.590 |
| 197 | 6.38 | 0.535 |
| 208 | 6.62 | 0.546 |
| 238 | 6.81 | 0.600 |

Other A>2 require both `--nuclear-radius-fm` and `--nuclear-surface-fm`.
These two flags also override the listed values. Both must be positive.
The same optional arguments are available on `electra ehijing run`; no local
or SLURM positional arguments or output layouts change. The SLURM submit script
accepts optional `NUCLEAR_TARGET=on|off`, `NUCLEON_WIDTH_FM`, `NUCLEAR_RADIUS_FM`
and `NUCLEAR_SURFACE_FM` key=value parameters. Python forwards the
options only when supplied, preserving invocation of older executables.

The [eHIJING model](https://arxiv.org/html/2304.10779v2) still supplies the gluon
TMD, transfer proposal/veto kernel, Qs equation, HT/GHT splitting kernels and
high/low-scale shower. Forward-ray Gaussian thickness replaces the continuous
medium thickness. Identity and collision time are sampled jointly from each
nucleon's contribution along that ray. Radiation gluons sample the remaining
medium starting at their existing formation-variable displacement. There is
no additional hard-sphere exclusion length: primary exclusion is by identity.
Qs uses the existing self-consistency solver outside the tabulated domain;
the high-scale separation retains a minimum of mu^2. The interpolation upper
endpoint is now safe without reading beyond the table.

All sampled elastic collisions are applied in nuclear-rest-frame time order.
Repeated hits, including hits from different partons, update the same whole
on-shell remnant. The solver retains the sampled transverse transfer, chooses
the forward two-body solution connected to zero transfer and conserves both
four-momentum and masses at that collision. A hit nucleon is split into one
quark–diquark pair after all hits, using its fixed proton/neutron species and
its final recoil momentum. Repeated hits exchange existing colour lines;
medium radiation inserts gluons into an existing supporting remnant's string.
No additional baryon is created for each collision or emission. Lund
fragmentation and its vertex machinery are retained. The pair's displacement
is the last sampled collision point, not its final momentum extrapolated back
to the collision time.

A transfer without a physical forward solution rejects the entire attempt and
regenerates its hard event and nucleus. Transfers are not clipped or resampled
inside an accepted event. Final B/Q must also match the untouched target
identities; a disagreement rejects the attempt. Only accepted attempts are
written. Logs report recoil and balance rejection counts. A run exits with a
diagnostic after 10000 consecutive unsuccessful attempts, rather than looping
forever. Rejection changes the accepted-event distribution; monitor its rate
before interpreting production results. These changes have synthetic tests,
but no physics production sample has been validated yet.

Exact conservation here applies to elastic recoil and remnant splitting.
The existing low-scale radiation step still subtracts its gluon spatial
momentum and restores the parent on shell independently; its approximate
energy balance has not been redesigned. This implementation therefore does
not claim exact global event energy conservation, nuclear excitation dynamics,
or reconstructed high-scale shower formation positions.

## Nuclear output bookkeeping

For neutral-current DIS, the writer counts signed final baryon number and
charge, including antibaryons, resonances, charged mesons and final decay
leptons. Intermediate particles and the scattered DIS lepton are excluded.
The required spectator counts are `Np = Z - Q_final` and
`Nn = A - B_final - Np`, computed in thirds to retain fractional quantum numbers.

In explicit-target mode these counts must equal the actual identities with
zero hits. The writer emits those identities at their original sampled
positions with rest-frame nucleon momenta, boosted to the output frame and
freestreamed to its t=0 hypersurface. It never deletes or invents spectators to
repair a target mismatch. The OSCAR header counts all emitted rows. Metadata
and the history's output snapshot contain a `nuclear_target` ledger with
identity, PDG species, hit count, sampled position in GeV^-1, primary ID and
geometry parameters. `nuclear_bookkeeping` reports quantum numbers and
`spectator_policy="untouched_target_identities"`.

With `--nuclear-target off`, admissible required counts set the old writer's
newly sampled spectators. Non-admissible counts retain the legacy A-1 fallback,
labelled `admissible=false` and `spectator_policy="legacy_unmodified"`. This
fallback does not apply to explicit-target events. Charged-current DIS remains
outside the supported hadronic charge budget.

OSCAR particle rows append `prod_time_status` after `xsecfac`: `0` means no
production reference is available (`prod_time` is a zero placeholder), `1`
means an available reference, including a genuine zero or an origin-only
assignment, and `2` means not applicable for pre-existing spectator nucleons.
Spectators retain `prod_time = form_time = 0`. Non-finite available production
references are errors.

The checked-in SMASH template declares
`Optional_Quantities: ["ID", "charge", "prod_time", "form_time", "xsecfac", "prod_time_status"]`.
The reader forces the profile start time to zero for status `0` and prints a
warning to stderr with the particle's PDG code; status `1` and `2` use the
supplied times without that warning. Warnings occur during startup validation
as well as when particles are imported. This explicitly retains the zero
fallback without assigning a physical production event to missing information.
Rebuild eHIJING and SMASH (or their container) to use the new interface. For
older OSCAR files, omit `prod_time_status` from `Optional_Quantities`; the reader
then retains the numerical-time behavior without availability information.
The local Python and SLURM workflows both use this template unless overridden;
custom configurations must match the columns and their order.

SMASH uses the imported Lund-based times by default (`Formation_Time.Mode:
OSCAR` under `Collision_Term.Cross_Section_Formation_Scaling`). The bundled
SMASH config omits this section to use the default implicitly. An explicit
secondary `Constant_Proper_Time` mode is available by setting `Mode` in that
section and specifying `Proper_Time` (default 1.0 fm/c).
It retains the OSCAR production time and replaces the endpoint by
`t_prod + (E/m) * Proper_Time`, with the proper interval in fm/c and status-2
spectators already formed. There is no automatic fallback from missing Lund
formation times to this mode.

## Diagnostic spacetime history

When enabled, each successful eHIJING event writes `event_00000042.history` next to its
OSCAR and metadata files, using the same eight-digit event ID and shard directory.
The file is native **HepMC3 Asciiv3 ASCII**, written with `HepMC3::WriterAscii`,
with momentum in GeV and vertex coordinates in mm (time in mm/c). It covers
eHIJING through output; it does not include SMASH transport. Collection defaults
to off. Enable it with `--history on` on either `ehijing` or `electra ehijing run`.
Existing local and SLURM callers use this default with the rebuilt executable.
The Python wrapper only passes this optional flag when explicitly requested,
so its default invocation remains compatible with older executables.

One file contains several HepMC event blocks representing diagnostic snapshots
of **one physical event**. `E` numbers identify snapshot sequence, starting at
zero. The `event_id` attribute preserves the physical event ID, including 64-bit
IDs; `attempt` identifies the successful generation attempt. Failed attempts
are discarded along with their event, rather than being attached to a later
successful event's history. Snapshot sequence is execution order, not physical
time. In particular, the high-scale shower is scale ordered and the existing
collision sampler retains its native order in attributes. Explicit elastic
recoil updates are subsequently sorted by their sampled rest-frame time.

The stages are:

- `hard_process`: the PYTHIA process record retained after generation.
- `nuclear_target`: the frozen sampled target, with `nuclear_target_id` on each
  nucleon and the primary identity in the comment.
- `medium_exact_recoil`: outgoing jet and whole remnant after each exact
  elastic collision; source shower index and target identity are separate.
- `high_scale_shower`: the complete PYTHIA event record before `Modified_FF`,
  including incoming particles, intermediate partons, remnants and leptons.
- `medium_parent_after_emission` and `medium_parent_after_kick`: observed parent
  momentum updates, with formation-variable or collision-index annotations.
- `medium_emission_products` and `medium_collision_products`: sampled gluons
  and recoil quark/diquark products, before append loses their parent fields.
  Indices identify their eventual entries in the shower record. Gluon collision
  samples are retained here even though the operative append does not retain them.
- `medium_appended` and `medium_modified`: newly assigned displacements and the
  complete record after medium processing.
- `hadronizer_import` and `hadronizer_input`: copied relative vertices and then
  absolute vertices, with explicit mappings back to source shower indices.
- `hadronization_and_decays`: the full hadronizer record, including intermediate
  hadrons and decayed particles, not just its final-state return vector.
- `oscar_output`: final hadrons extrapolated to t = 0 and the actual spectator
  nucleons retained from the explicit target (or sampled by the legacy writer). Coordinates are recorded before the
  OSCAR file rounds them to five decimals in fm.

Only stages with observed medium updates/products are emitted. Identity is
`(record, record_index)`; indices are local to `process`, `shower`, `hadronizer`
or `oscar`. `source_record`/`source_index` link imports and appended products;
`source_index=-1` means no direct import association. `nuclear_target_id`
identifies a target nucleon for nuclear observations, recoil pairs and imported
primary remnants; `-1` means unassigned. Shower snapshots also retain the
separation scale and primary target identity. Native `mother1`,
`mother2`, `daughter1`, `daughter2`, `pythia_status`, colour tags, scale and
sampled collision arrays are attributes. PYTHIA's range and copy conventions
still apply to those ancestry fields. Final hadronizer particles have an
`oscar_id`; OSCAR observations map back to the corresponding hadronizer index.
Spectators are marked separately and do not acquire invented shower ancestry.
The output snapshot also retains the actual OSCAR production/formation times,
`oscar_prod_time_status`, and cross-section factors, alongside its extrapolated
t = 0 positions. Native
decay endpoints from PYTHIA's stored vertex and lifetime are recorded when
available in the hadronizer record.
The PYTHIA system entry at index zero is retained and marked as bookkeeping.

Each simulated particle has its own observation vertex, with status 900.
These are **observations, not newly inferred physical interaction vertices**.
To retain compatibility with HepMC3 readers that reject vertices with empty
incoming lists, each has one auxiliary incoming particle (PDG 0, status 4,
`diagnostic_anchor=1`). Filter these serialization anchors when counting
simulated particles. They carry no physics meaning. Physical ancestry is in
the attributes rather than invented connections between observation vertices.
Standard event displays will therefore show an observation graph, not a
reconstructed shower worldline graph.

`has_vertex=0` means a stored zero is unavailable; `has_vertex=1` distinguishes
an explicitly assigned zero. `stored_vertex_xyzt` retains the raw value with
`stored_vertex_units`, while HepMC `V` positions use mm. All four-coordinate
attributes use `(x,y,z,ct)`. Nuclear geometry origins are separately recorded
in GeV^-1, along with the target momentum defining the nuclear rest frame.
The collector adds no boosts or translations. Explicit nuclear/recoil stages
already supply boosted positions and momenta from the model; other snapshots
retain their stage's existing frame and origin conventions. The supplied workflow
requires `PartonVertex:setVertex=off`; custom settings enabling it still require
the import-policy review described in [Unit conventions](docs/units.md).

Production references and geometric breakup pairs are retained as attributes
when present. The native fragmentation trace (including Lund geometry and
special string branches) is embedded in `comment_fragmentation`, even when the
separate `.fragmentation` file is disabled. Comments are HepMC string attributes
(`A` records), with newlines escaped by the ASCII writer, so standard readers
can read the file without stripping custom text lines.

Collection reads/copies existing values and consumes no random draws. It does
not restore operative ancestry, assign shower lifetimes, move vertices, or
change momenta, colours, decays, collision sampling or OSCAR output. Missing
spacetime information remains visible as a gap. File size and memory usage grow
with the full particle records and observed medium updates; leave collection off
for production runs that do not need these diagnostics.

After building with `BUILD_TESTING=ON`, run `ctest --test-dir <build-dir>
--output-on-failure`. The synthetic history test checks native HepMC3 reader
round trips, identity mappings, collision ordering, missing versus known-zero
vertices, units, multiline comments and identical OSCAR/metadata output with
diagnostics enabled. It does not generate simulation events.

## Hadronization spacetime trace

Enable diagnostics in the hadronization settings file passed to
`--hadronization-config`:

```text
HadronVertex:trace = on
Fragmentation:setVertices = on
HadronVertex:mode = 2
```

`HadronVertex:trace` defaults to off. When enabled, each successful event gets
an `event_XXXXXXXX.fragmentation` file beside its `event_XXXXXXXX.oscar` and
metadata files, in the same shard directory. The file contains only the successful
generation attempt, including its internal hadronization passes. Failed generation
attempts are discarded. Error messages still go to `stderr`.

Rebuild the bundled PYTHIA library and the eHIJING executable and use the updated
bundled XML settings. A container/SIF must include all of these updates.

The dump identifies event, generation attempt, hadronization pass and string
subsystem. It prints ordered partons and color connections, local breakup
vertices, inherited origins, absolute breakup pairs, prescription displacements,
saved vertices and the final complete event record. Coordinates are ordered
`(t,x,y,z)`, with time in fm/c, position in fm, and momentum in GeV. Section
headings, particle names and an exported-hadron time summary make the text easier
to scan; values use 12 significant digits with scientific notation when needed.
In mode 2 the saved vertex is the formation point. Other modes retain their
existing interpretation. Version 3 adds fragmentation-side labels, momenta at
fragmentation, and `GEOMETRY`, `REGION`, and `PARAMETER_BREAKUP` records. These
describe PYTHIA's longitudinal massless-reference parameter plane; diagnostics
fill unused regions in a copy of the system without changing the event or RNG.
Version 4 additionally records the explicit inherited origin and
`REFERENCE_BREAKUP` points after spacetime conversion corrections, before heavy
endpoint offsets and vertex smearing. These include special finalTwo solutions.
The trace preserves raw fractions alongside the converted points.

Final entries include the OSCAR `output_ID` for exported particles, the source
string and branch, and indicators showing whether their saved vertex or breakup
metadata changed after fragmentation. Hadrons created later by decays or copies
list ancestor strings through their parent links; multiple sources are retained.
One-hadron collapses use their inherited origin as the production point without
constructing a breakup pair. Full production-reference vertices are recorded
separately from geometric breakups/turns. The
dump precedes the writer's freestreaming to t = 0. The writer adds no spatial
smearing to hadron positions.
Diagnostics do not change particle data or the OSCAR columns.

### Plotting a trace

ELECTRA provides `utils/plot_fragmentation.py` (Python 3.10+ and matplotlib).
From the ELECTRA repository root:

```bash
python3 utils/plot_fragmentation.py /path/to/event_00000000.fragmentation
python3 utils/plot_fragmentation.py /path/to/event_00000000.fragmentation \
    --string 0 --pass 1 --axis z -o string_0.pdf
python3 utils/plot_fragmentation.py /path/to/event_00000000.fragmentation \
    --string 0 --pass 1 --frame string-rest -o string_0_rest.pdf
```

The default PDF sits beside the input file and contains the common lab-frame
overview, individual string/pass diagrams, Lund reference pages, and
final particle-time tables in the lab and each available input-string rest
frame by default (`--no-table` omits them). The container script
is at `/opt/electra/share/utils/plot_fragmentation.py`.

Lund boundaries are enabled by default on the overview and saved-vertex pages,
including `--overview-only`. Solid bold lines show original endpoint boundaries,
thin lines show breakup constituents, and verified reference string area is
shaded. Filled green and purple circles show stored breakup and formation
vertices. Boundary turns remain in the curves without separate markers or
labels. Hollow `RB`/`RF` circles show converted reference breakups/meetings.
The primary vertex is a yellow star. Hadron worldlines use `dx/dt = p/E` from
the saved formation point; they exclude later
interactions. Nonfinal particles are not assigned an invented lifetime.

Reference paths use the explicit inherited string/leg origin, never an average
of endpoint birth events. Canonical half-gluon momenta normalize rescaled
regional bases before stitching. Spacetime conversion corrections are included;
paths and area fill must pass domain, ordering and continuity checks. `B!` marks
a converted point outside the classical reference domain; `H!` marks a stored
lab-frame formation before an associated breakup. Corrections in PYTHIA can
produce such assignments; no causal classical connection is fabricated.

Reference meeting points exclude transverse fragmentation momentum, while the
saved late formula includes the full hadron four-momentum. Therefore turning
off `HadronVertex:smearOn` does not make these points identical in the lab. The
reference picture also excludes heavy-endpoint offsets and transverse tunnelling.
See [the late prescription, Eq. (25)](https://arxiv.org/html/1808.04619v2#S3.SS1).

### Single-hadron collapses

The default `--collapse-model both` overlays a labelled one-hadron Lund
diamond and the original massive input endpoint curves, also on the dedicated
collapse page. Dash-dot lines are the effective diamond; dotted lines are the
input curves. The diamond uses the final hadron momentum after recoil exchange
and closes at
the stored late point `F = O + p_h/kappa` and uses the positive input endpoint
direction in the final hadron rest frame as its axis convention. The proper
model extent is `m_h/kappa`. Turns remain in the curves without markers or labels.
This is a massless effective post-collapse string, not the physical evolution
of the original massive partons. It is used when the assigned vertex satisfies
that formula; smearing/heavy corrections can require the input-model fallback.

`--collapse-model effective` selects only the diamond; `--collapse-model input`
selects only the massive constant-force diagnostic. Each
endpoint starts at its own recorded birth vertex, turns at
`abs(p_parallel)/kappa` and stops on its first return at
`2*abs(p_parallel)/kappa`, in the input string CM frame. It excludes pair creation
and recoil exchange, so it need not meet the saved hadron formation. See
[the massive equation of motion, Eq. (17)](https://arxiv.org/html/1808.04619v2#S2.SS2.SSS5).

Recorded starts (hollow diamonds) describe input constituent availability, not
hadron formation. Neither input endpoint turns nor effective diamond turns create
constituents. The effective origin is an inherited model anchor, while the purple
point is the assigned hadron formation. No recoil-exchange history is reconstructed.

### Frames and compatibility

All main plots for an event and frame/spatial projection use shared axis limits.
Longitudinal pages share numerical ranges too, while retaining each string's
own CM and axis; their coordinates cannot be overlaid as one physical frame.
Endpoint-detail insets retain their local zoom. Reference projection pages
use the same boost and coordinate origin as the corresponding saved pages.

Use `--frame lab --overview-only` to plot all selected strings in one common
frame. `--string 0 --pass 1` restricts the selection. Reference pages default
to each string's own CM longitudinal axis; effective collapse pages use the
final hadron rest frame. `--lund-view projection` instead uses the chosen
`--frame`/`--axis` on reference pages. The default `--frame both` includes lab
and input-string CM saved-vertex plots and tables. A multi-string overview
uses only lab; each individual string has its own rest view when its momentum
is timelike. `--frame string-rest` requests those individual rest views and
keeps the multi-string overview in lab. `--frame lab` uses lab reference pages
too unless `--lund-view string-rest` is explicitly selected. `--overview-only`
omits per-string pages and tables.

Tables show `prod_time` and a single `form_time` in mode 2, using the final particle
vertex when available, otherwise the fragmentation assignment. They recompute
rest-frame production time from the full production references (breakups or
origins). Old missing production times remain `unavailable`, distinct from zero.
Rest tables use the input string CM and primary origin, matching saved-vertex
rest plots; effective collapse reference pages use the final hadron CM and
inherited origin.
Rest-table rows are grouped by string and collected on the same page; each row
uses its own string CM, not a common rest frame for the event.

Individual string rest frames cannot be mixed
into a common physical spacetime overview.

Other options include `--axis x`, `--axis y`, `--no-labels`, `-o output.pdf`,
and `--no-lund-boundaries` to omit the overview overlay. `--duration 5` extends
free hadron lines 5 fm/c past the latest formation; it does not set constituent
trajectory durations.

The plotter reads v2/v3/v4. Old unsmeared light-flavour breakup pairs provide
converted reference points when v4 records are absent. Rebuild bundled PYTHIA,
eHIJING and the container, then regenerate events to apply the native fixes:
collapse recoil reshuffling changes momentum while preserving past vertices,
and both hadrons of a two-hadron mini-string inherit the positive endpoint
origin, matching ordinary fragmentation. This origin is a consistent inherited
anchor, not a newly determined colour-exchange event. Original endpoint births
remain in the trace. Old traces preserve their historical metadata, including
any recoil-boosted vertices or different mini-string inherited origins.

### Production references and collapse origins

For ordinary open strings in mode 2, production references for the first and
last hadrons replace the corresponding endpoint turn with the inherited string
origin. Interior hadrons retain both breakup references. Geometric boundary
points remain available for the unchanged late formation prescription and Lund
paths. Closed strings, junctions and two-hadron mini-strings retain their previous
production rules. Production references and any inherited origin transform with
the particle under genuine rotations/boosts and survive particle copies.

A single-hadron collapse has an origin-only production assignment. The existing
native origin is

\[
O^\mu=\tfrac12\left(v_{\rm first}^\mu+v_{\rm last}^\mu\right)
\]

when the first endpoint has a stored vertex; otherwise the default origin is
zero. This is an inherited average, not a measured colour-exchange vertex.
The native formation assignment is

\[
F^\mu=O^\mu+S^\mu+r\,p_h^\mu/\kappa,
\qquad t_{\rm prod}=O^0,
\]

where `S` is the configured smearing displacement and `r` is the native heavy
flavour `redOsc` factor. For unsmeared light-flavour collapses, `S=0`, `r=1`,
so `F=O+p_h/kappa` and `t_form=t_prod+E_h/kappa`. The momentum is the final
hadron momentum after recoil exchange. These equations use coordinates in fm
and momentum in GeV, with `kappa` in GeV/fm; stored PYTHIA vertices are converted
onsite to mm. The production update does not modify the formation calculation.
Rebuild PYTHIA/eHIJING and regenerate events to apply the new production rules;
historical traces retain their recorded times.

## Reproducibility and seeds

`--seed` is an ELECTRA master seed: an unsigned decimal integer in
`[0, 18446744073709551615]`. It is not a raw PYTHIA seed. In particular, master
seed zero is reproducible. The Python wrapper passes it unchanged; omitting it
in standalone eHIJING obtains entropy once and prints the selected master seed.
Rebuild both the bundled modified PYTHIA library (including its XML settings)
and eHIJING to use this policy. A linked policy-version symbol prevents silently
using an older medium library. Previously generated samples are not reproduced
by the new policy; retain their original binaries and records.

The guarantee is complete event replay in fresh processes for the same master,
chunk starting ID, settings, table bytes, executable and runtime environment.
Streams advance through failed attempts and accepted events without reseeding.
The recorded successful attempt number counts all attempts for that output
event, including across the driver's 15-attempt batches. Changing the requested
event count or output shard width preserves the event prefix for an unchanged
chunk start. Changing chunk starts changes the sample. PYTHIA can retain and
adapt sampling state, so reseeding by event ID would not suffice for independent
event replay. Diagnostic collection consumes no additional draws.

The frozen `electra-seeds-v1` derivation is implemented in
`external/pythia8/include/eHIJING/SeedPolicy.h`. Define `mix(x)` by unsigned
64-bit arithmetic:

```text
x += 0x9e3779b97f4a7c15
x = (x XOR (x >> 30)) * 0xbf58476d1ce4e5b9
x = (x XOR (x >> 27)) * 0x94d049bb133111eb
return x XOR (x >> 31)
context = mix(master XOR mix(first_event_id) XOR 0x454c454354524131)
stream_seed = mix(context XOR mix(stream_id))
```

Each `mt19937` receives `seed_seq{low32(stream_seed), high32(stream_seed)}`.
Each PYTHIA instance receives `1 + stream_seed % 900000000`, with
`Random:setSeed = on`, after reading settings and before `init()`. PYTHIA's own
seed zero selects the time and negative seeds select its default sequence;
neither reaches PYTHIA through the ELECTRA master-seed conversion.

| ID | Stream owner | Use |
| --- | --- | --- |
| 1 | Hard PYTHIA | Hard process, ordinary shower draws, explicit nuclear target |
| 2 | PartonLevel geometry | Legacy hard interaction position |
| 3 | TimeShower geometry | Deterministic path/thickness calculations; no draws currently |
| 4 | TimeShower medium, MultipleCollision base | High-scale collisions and target selection |
| 5 | TimeShower medium, eHIJING derived class | High-scale radiation |
| 6 | Modified_FF | Low-scale emissions and remnant sampling |
| 7 | Modified_FF collision sampler | Low-scale additional collisions |
| 8 | Modified_FF geometry | Deterministic calculations; no draws currently |
| 9 | Modified_FF colour-choice engine | Equal-probability colour/recoil choices formerly using C `rand()` |
| 10 | Hadronizer PYTHIA | Fragmentation and decays |
| 11 | Hadronizer helper engine | Legacy nucleon and remnant choices |
| 12 | Event writer | Legacy spectator positions; unused for explicit targets |
| 13–16 | Hadronizer PartonLevel/TimeShower medium objects | Separate namespace; inactive in the supplied hadronization workflow |

The unused driver engine was removed. The unused `Event` default engine and
uncalled GSL VEGAS helper are not operative streams. The fixed stream labels do
not guarantee disjoint sequences: mapping into PYTHIA's finite seed range can
produce collisions. Cross-platform bitwise equality and matched random histories
between changed physics settings are outside this contract.

Per-event metadata's `reproducibility` object records policy, master, chunk start,
successful attempt, effective PYTHIA seeds and the 64-bit inputs to `seed_seq` for
other engines. Wide seed values are decimal strings to retain all bits in JSON
readers. `inputs_sha256` identifies settings and Qs/GHT file bytes. Histories
include the same provenance as an output-observation attribute. The fragmentation
configuration label uses the canonical settings path, so alternate spellings of
the same file agree; relocating that file changes its diagnostic label.

New tables use `max_digits10` so in-memory values match reloads. Generation is
locked per table and publishes complete files atomically. Existing cached tables
are loaded without rewriting or requiring cache-directory write access. Retain
the same table bytes for replay; regenerating historical low-precision tables can
change numerical inputs.

The CTest `seed_policy` test runs without event generation. The parent ELECTRA
repository also provides an explicit smoke check with small temporary runs:

```bash
PYTHONPATH=python python3 tests/reproducibility_smoke.py --executable /path/to/ehijing
```

It compares complete records in fresh processes for HT/GHT and explicit/legacy
targets, generated/reloaded tables, history on/off, direct/wrapper runs and
prefix preservation. Its outputs are temporary; research/reference data and
installed utilities are untouched.

## Generating documentation

Install Doxygen first. Installing Graphviz is also recommended because the
configuration enables class diagrams through `dot`.

Generate HTML documentation from this directory with:

```bash
doxygen Doxyfile
```

or through CMake:

```bash
cmake -S . -B build
cmake --build build --target docs
```

To configure the documentation target without requiring the full eHIJING
physics dependencies, use:

```bash
cmake -S . -B build-docs -DEHIJING_BUILD_DOCS_ONLY=ON
cmake --build build-docs --target docs
```

The generated documentation is written to `docs/doxygen/html/index.html`.
