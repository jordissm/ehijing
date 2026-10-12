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
In mode 2 the saved vertex is the formation
point. Other modes retain their existing interpretation.

Final entries include the OSCAR `output_ID` for exported particles, the source
string and branch, and indicators showing whether their saved vertex or breakup
metadata changed after fragmentation. Hadrons created later by decays or copies
list ancestor strings through their parent links; multiple sources are retained.
One-hadron collapses explicitly report production time as unavailable because
no breakup pair is reconstructed. The
dump precedes the writer's freestreaming to t = 0 and output spatial smearing.
It describes computed vertices and their connections, rather than reconstructing
continuous string trajectories. Diagnostics do not change particle data or the
OSCAR columns.

### Plotting a trace

ELECTRA provides `utils/plot_fragmentation.py` (Python 3.10+ and matplotlib).
From the ELECTRA repository root:

```bash
python3 utils/plot_fragmentation.py /path/to/event_00000000.fragmentation
python3 utils/plot_fragmentation.py /path/to/event_00000000.fragmentation \
    --string 0 --pass 1 --axis z -o string_0.pdf
```

The default PDF sits beside the input file and includes a time-versus-z overview,
individual string/pass diagrams and hadron time tables. Use `--axis x` or
`--axis y` for other projections, `--no-labels` for crowded diagrams, or
`--overview-only` for a single page. In the ELECTRA container the script is at
`/opt/electra/share/utils/plot_fragmentation.py`.

Plots use absolute breakup and saved coordinates before OSCAR freestreaming and
transverse smearing. Lines connecting breakup and saved vertices are geometric
guides, not particle worldlines. An earlier breakup is shown even when its
production time was not stored; the time table marks missing values as
`unavailable`. Faded hadrons were not exported, and changed final vertices have
separate square markers. Later decay/copy hadrons use their ancestry when it
identifies a single string; ambiguous ancestry is shown as unassigned.

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
