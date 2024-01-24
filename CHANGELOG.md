# Changelog

## NEXT

- Create `install.sh` to replace `setup.sh` with some features moved from `build.sh`
- Add reference radius and reference energy to runtime options
- Allow user to pass ideal shock angles (i.e., theta, phi, and width) in degrees
- Allow user to pass observer angles (i.e., theta and phi) in degrees
- Change value of `VERSION` to be consistent with this and configure files
- Refactor parameter-related constants
- Improve config-file read and echo
- Always output the initial time step
- Accept different widths for the ideal shock
- Define a merged shock-scale parameter

## v0.3.0 (18Dec2023)

- Implement point-observer output
- Refactor and add features to build-process tools

## v0.2.6 (03Aug2023)

- Add `--download-ext-deps` to `setup.sh` CLI.

## v0.2.5 (11Apr2023)

- Redefined options for `setup.sh`.

## v0.2.4 (06Jan2023)

- Allow `idealShockFalloff` to be 0 and change default value to 0.

## v0.2.3 (20Dec2022)

- Expanded installation instructions in README.
- Added `--with-ext-deps=DIR` to configure options.
- Moved `--with`-style configure options to `setup.sh`.
- Copy `config.log` during `make install`.

## v0.2.2 (11Nov2022)

- Revert to multiplying velocity by the ideal shock factor because EPREM nodes are not in the frame co-moving with the shock.

## v0.2.1 (11Nov2022)

- Avoid running `autoreconf` as part of setup.sh if possible.

## v0.2.0 (10Nov2022)

- Exponentially relax shocked quantities to unshocked values downstream of shock.
- Allow user to control rate of relaxation with `idealShockFalloff`.

## v0.1.1 (09Nov2022)

- Fix velocity scaling for ideal shock in `flow.c`.

## v0.1.0 (02Nov2022)

- Initial release of uncoupled EPREM
