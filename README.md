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

Final entries include the OSCAR `output_ID` for exported particles, the source
string and branch, and indicators showing whether their saved vertex or breakup
metadata changed after fragmentation. Hadrons created later by decays or copies
list ancestor strings through their parent links; multiple sources are retained.
One-hadron collapses explicitly report production time as unavailable because
no breakup pair is reconstructed. The
dump precedes the writer's freestreaming to t = 0 and output spatial smearing.
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

The default PDF sits beside the input file and includes a time-versus-z overview,
individual string/pass diagrams, longitudinal constituent worldline pages for
v3 geometry records, and hadron time tables. Use `--axis x` or
`--axis y` for other projections, `--no-labels` for crowded diagrams, or
`--overview-only` for a single page. `--frame string-rest` boosts both vertices
and momenta to one selected string's rest frame, with the primary vertex as the
coordinate origin. `--duration 5` extends outgoing hadron worldlines to 5 fm/c
after the latest formation point in the plot frame. In the ELECTRA container the script is at
`/opt/electra/share/utils/plot_fragmentation.py`.

The Lund pages default to each string's own rest frame, with `z_parallel`
aligned with the recorded positive light-cone direction. A two-hadron mini-string
may choose this reference decay axis differently from its input parton axis.
`--lund-view projection` instead
uses the selected `--frame` and `--axis` on those pages. This option does not
change the saved-vertex plots or time tables.

Breakup vertices are solid green circles, formation vertices solid purple
circles, and the primary vertex a yellow star. Saved-vertex diagrams show
free outgoing hadron worldlines with `dx/dt = p/E`, beginning at the stored
formation point; they do not connect adjacent breakup points to that vertex.
These lines are extrapolations without downstream interactions or decays.

Saved-vertex pages also include massive original endpoint curves for two-parton
open strings, including one-hadron collapses. The dash-dot curves integrate
`dp_parallel/dt = +/- kappa` and `dx/dt = p/E` in the input string CM frame,
using the recorded on-shell masses and momenta, and boost the spacetime paths
back to the selected frame. Diamonds mark recorded endpoint starts; triangles
mark model turns. Curves show one outward-and-return excursion and stop there.
Individual string pages include an endpoint close-up when the motion is small
compared with the formation-time scale, and diagram labels are spaced apart.
This is an ideal unbroken constant-tension model, not an exact constituent
history from PYTHIA. Each endpoint starts at its own recorded vertex with its
own force clock; distinct birth vertices need not describe a single coherent
straight string. Fragmentation/pair creation and collapse recoil exchange are
not modelled, and the curves are not fitted to saved formation points. Massive
paths through gluon kinks/junctions are omitted with a diagnostic message.
The model uses [the straight-string equation of motion, Eq. (17)](https://arxiv.org/html/1808.04619v2#S2.SS2.SSS5).
New traces explicitly record each input parton's `mass_GeV`; older v2/v3 files
can recover that mass from the recorded four-momentum.

The separate longitudinal pages follow the classical massless-reference Lund
constituent characteristics, splitting paths at gluon-region boundaries and
including outward endpoint motion before the boundary turns. Original
endpoint boundaries are bold solid lines, breakup-constituent boundaries are
thin solid lines, and the surviving string worldsheet is shaded. Pair creation
removes the future quadrant of each breakup in the string parameter plane;
the shading excludes those regions. If any pair-creation coordinates are
missing, the area is not filled. Outward endpoint arms are reconstructed
independently of hadron solutions so missing final-two coordinates do not
hide an otherwise recorded original endpoint arm. Their green and
purple points, and outgoing longitudinal hadron slopes, belong to that reference
picture. They exclude transverse tunnelling/smearing and heavy-endpoint
corrections and can differ from the saved physical vertices. PYTHIA's smeared
late prescription need not lie in both breakup vertices' future light cones;
straight connections would therefore not give physical constituent worldlines.
Missing final-two coordinates or discontinuous/invalid regions are reported
instead of interpolated. A shared mean endpoint vertex (or mean origin for a
leg without recorded endpoints) anchors all reference paths of a string/leg,
so neighbouring paths join; saved hadron metadata is unchanged. One-hadron
collapses have no parameter-plane paths for newly created constituents, but
their original massive endpoints are shown by the model above.

Old v2 files plot saved vertices, massive input endpoints and hadron worldlines;
regenerate traces after rebuilding PYTHIA/eHIJING to obtain the parameter-plane
constituent pages. Time tables mark
missing production times as `unavailable`. Faded hadrons were not exported,
and changed final vertices have square markers. Later decay/copy hadrons use
their ancestry when it identifies a single string; ambiguous ancestry is
shown as unassigned. Nonfinal hadrons are not extrapolated with an invented lifetime.

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
