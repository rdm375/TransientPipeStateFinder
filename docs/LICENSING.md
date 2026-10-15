# Licensing and dependency provenance

## TransientPipeStateFinder

TransientPipeStateFinder is licensed under the GNU Affero General
Public License, version 3 only (AGPL-3.0-only).

The complete license text is provided in the repository's LICENSE file.

## Forward solver

TransientPipeStateFinder uses the independently developed
isothermal transient-pipe forward solver:

Repository:
https://github.com/rdm375/transient-pipe-flow

License:
AGPL-3.0-only

Licensed main-branch revision:
8eba6e7

Previously validated numerical baseline:
670ac493e4461e678ca66f2254dc3d402df3f312

The licensing commits added documentation and the license text,
without modifying the numerical implementation.

Historical validation results retain their original source revision.
Future validation and engineering reports must identify the actual
solver revision used to build the executable.

## Third-party dependencies

The Git history identifies Richard Myers as the recorded author of
the forward solver's commits, and no Git submodules were identified.

These observations do not independently establish the licensing
status of every external dependency or copied component.

Third-party code, libraries, and distribution obligations must be
reviewed before distributing combined binaries.

## Reproducibility

Engineering reports should record:

- StateFinder source revision and modification status.
- Forward solver source revision and modification status.
- Compiler and build configuration.
- Numerical discretization and solver settings.
- Measurement configuration and uncertainty assumptions.
- Forecast horizon and engineering acceptance tolerances.

Source revisions must correspond to the code actually compiled,
rather than the repository state when a report is generated.
