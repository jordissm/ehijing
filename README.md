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

`HadronVertex:trace` defaults to off. The trace is written to `stderr`; redirect
it with `2> hadronization.log` on the eHIJING command. Other error messages also
appear in that file. Rebuild the bundled PYTHIA library and the eHIJING executable
and use the updated bundled XML settings before enabling this new option.

The dump identifies event, generation attempt, hadronization pass and string
subsystem. It prints ordered partons and color connections, local breakup
vertices, inherited origins, absolute breakup pairs, prescription displacements,
saved vertices and the final complete event record. Coordinates are ordered
`(t,x,y,z)`, with time in fm/c, position in fm, and momentum in GeV; scientific
notation preserves small values. In mode 2 the saved vertex is the formation
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
