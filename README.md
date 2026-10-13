## Documentation

The eHIJING headers and sources are documented with Doxygen.

[Unit conventions](docs/units.md) describe spacetime units, conversion factors,
interface requirements, and remaining review items.

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
