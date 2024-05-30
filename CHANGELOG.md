# Changelog

## Developer Notes

When incrementing the version number to X.Y.Z, please do the following
* create a new subsection here (below **NEXT**) with the title vX.Y.Z (YYYY-MM-DD)
* update the version number in `configure.ac` and in `src/global.h`
* commit with the message "Increment version to X.Y.Z"
* create a tag named "vX.Y.Z" with the message "version X.Y.Z"
* push and follow tags

## NEXT

## v0.5.0 (2024-05-30)

- Add command-line option to print version number and exit
- Replace `obsUseDegrees` and `idealShockUseDegrees` with single `useDegrees` parameter
- Add `DEG2RAD` and `RAD2DEG` constants
- Fix bug in `make dist`: add missing `tools.sh`
- Remove redundant parameter descriptions in `install.sh` help text
- Fix bug in default `numSpecies`, `mass`, and `charge`
- Fix error in how `mhdB` and `parkerB` compute ideal-shock magnetic field

## v0.4.0 (2024-01-25)

- Create `install.sh` to replace `setup.sh` with some features moved from `build.sh`
- Change default value of `boundaryFunctBeta` from 1.7 to 2.0
- Add reference radius and reference energy to runtime options
- Allow user to pass ideal shock angles (i.e., theta, phi, and width) in degrees
- Allow user to pass observer angles (i.e., theta and phi) in degrees
- Change value of `VERSION` to be consistent with this and configure files
- Refactor parameter-related constants
- Improve config-file read and echo
- Always output the initial time step
- Accept different widths for the ideal shock
- Define a merged shock-scale parameter
- Fix redundant installation of external dependencies

## v0.3.0 (2023-12-18)

- Implement major algorithmic updates to adiabatic change and focusing developed by PSI
- Remove species dependence of "grid" arrays for energy, velocity, and momentum
- Set the default value of `rScale` to the value of 1 solar radius in au (i.e., `RSAU` as defined in `global.h`).
- Set the background (a.k.a "seed") spectrum minimum value to `DBL_MIN`
- Implement point-observer output
- Add initial components for model-agnostic MHD coupling
- Refactor and add features to build-process tools

## v0.2.6 (2023-08-03)

- Add `--download-ext-deps` to `setup.sh` CLI.

## v0.2.5 (2023-04-11)

- Redefined options for `setup.sh`.

## v0.2.4 (2023-01-06)

- Allow `idealShockFalloff` to be 0 and change default value to 0.

## v0.2.3 (2022-12-20)

- Expanded installation instructions in README.
- Added `--with-ext-deps=DIR` to configure options.
- Moved `--with`-style configure options to `setup.sh`.
- Copy `config.log` during `make install`.

## v0.2.2 (2022-11-11)

- Revert to multiplying velocity by the ideal shock factor because EPREM nodes are not in the frame co-moving with the shock.

## v0.2.1 (2022-11-11)

- Avoid running `autoreconf` as part of setup.sh if possible.

## v0.2.0 (2022-11-10)

- Exponentially relax shocked quantities to unshocked values downstream of shock.
- Allow user to control rate of relaxation with `idealShockFalloff`.

## v0.1.1 (2022-11-09)

- Fix velocity scaling for ideal shock in `flow.c`.

## v0.1.0 (2022-11-02)

- Initial release of uncoupled EPREM
