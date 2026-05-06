## Documentation

The eHIJING headers and sources are documented with Doxygen.

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
