# eHIJING Documentation

eHIJING is the event-generation component that combines PYTHIA DIS event
generation with eHIJING medium modification, hadronization, and OSCAR2013-style
output.

The executable workflow is:

1. Parse runtime configuration into `RunConfig`.
2. Configure and initialize the PYTHIA hard process.
3. Reconstruct `DISKinematics` and apply runtime `DISCuts`.
4. Apply low-\f$ Q^2 \f$ medium modification with `Modified_FF`.
5. Convert the modified partonic event into hadrons with `Hadronizer`.
6. Write OSCAR2013 event records and JSON metadata.

## Runtime Inputs

The generator is configured through command-line arguments and input files
rather than recompilation-sensitive constants.

- Hard-process PYTHIA settings are loaded from `RunConfig::hard_process_config_path`.
- Hadronization PYTHIA settings are loaded from `RunConfig::hadronization_config_path`.
- DIS event-selection cuts are loaded from `RunConfig::dis_cuts_config_path`.
- Tabulated eHIJING data are read from `RunConfig::tabulation_path`.

The DIS cuts file is parsed by `load_dis_cuts()` and expects the keys `yMin`,
`yMax`, `xBMin`, `xBMax`, `Q2Min`, and `W2Min`.

## Public API Map

- `RunConfig` and `parse_args()` define the executable runtime interface.
- `DISKinematics`, `DISCuts`, `compute_dis_kinematics()`, and `trigger()` define
  DIS reconstruction and event selection.
- `Modified_FF` owns the low-\f$ Q^2 \f$ medium-modification sampler.
- `Hadronizer` owns the dedicated PYTHIA instance used for string fragmentation.
- `write_event_headers()` and `write_event_output()` write OSCAR2013 particle
  records and per-event JSON metadata.
- `make_event_paths()` and related helpers define deterministic sharded output
  paths for large productions.
- `ehijing::constants` centralizes compile-time conventions that are not meant
  to change run by run.

## Conventions

The code relies on a few explicit PYTHIA/eHIJING conventions:

- PYTHIA DIS event-record indices are defined in `ehijing::constants::pythia`.
- PDG Monte Carlo particle IDs live in `ehijing::constants::pdg` and use the
  `_id` suffix.
- Masses and numerical cutoffs that are still compile-time constants are grouped
  in `ehijing::constants::mass`, `ehijing::constants::numeric`, and related
  namespaces.

## Generating The HTML

From `src/ehijing`, either run Doxygen directly:

```bash
doxygen Doxyfile
```

or use the CMake target:

```bash
cmake -S . -B build
cmake --build build --target docs
```

If you only want the documentation target and do not need to configure the
generator executable, use:

```bash
cmake -S . -B build-docs -DEHIJING_BUILD_DOCS_ONLY=ON
cmake --build build-docs --target docs
```

The generated entry point is:

```text
docs/doxygen/html/index.html
```
