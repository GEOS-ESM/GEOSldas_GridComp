# Changelog for `GEOSldas_GridComp`

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

-----------------------------

## [Unreleased]

### Added

### Changed

### Fixed

### Removed

### Deprecated

-----------------------------

## [v4.0.0] - 2026-09-10

- Generally not 0-diff vs. v3.2.0 (owing to revised QC of Tb, sfds, sfmc; also requires newer, non-0-diff GEOSgcm_GridComp).

### Added

- Added support for river routing, incl. ensemble simulations ([PR #145](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/145), [PR #174](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/174), [PR #176](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/176)).
- Added support for lake tiles (single ensemble member only) ([PR #181](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/181), [PR #176](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/176)).
- Added support for running ISSM (Ice-Sheet and Sea-level System Model; single ensemble member only) ([PR #161](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/161), [PR #176](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/176)). 

- Added assimilation of surface soil moisture observations from H-SAF ASCAT H121 CDR v8 and H139 ICDR netcdf products (MetOp-A/B/C) ([PR #186](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/186)).

- Added optional NetCDF4 output of ObsFcstAna; changed namelist variable "out_ObsFcstAna" from logical to integer ([PR #163](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/163), [PR #185](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/185)).

- Added Matlab and python readers for binary Tb scaling parameters files ([PR #179](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/179), [PR #191](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/191)).
- Added python reader for binary catparam files ([PR #191](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/191)).

- Added SMOS Tb preprocessing scripts ([PR #189](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/189)).


### Changed

- Added peatland QC for sfds and sfmc observations ([PR #186](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/186)).
- Added QC of SMAP L1C_TB using max value for Tb_error ([PR #190](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/190)).

- Renamed './cat' output directory to './diag'; created link from './cat' to './diag' for backward compatibility ([PR #198](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/198)).


- Revised and cleaned up RESTART options ([PR #160](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/160), [PR #166](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/166)):
  - Clarified scope and constraints of RESTART=1 and RESTART=2.
  - Added RESTART=3 (formerly RESTART=G, which had been removed).
  - Cleaned up RESTART=M.
- Updated Landice ("glc") HISTORY Collection to that of M21C plus key ISSM outputs ([PR #181](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/181)).
- Replaced legacy HDF4 Fortran interface with a C bridge and `ISO_C_BINDING` module ([PR #194](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/194)).
- Updated CI ([PR #181](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/181)).


### Fixed

- Fixed `read_obs_param()` parsing for the current obsparam format by reading forecast variable names and units ([PR #185](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/185)).
- Fixed crashes in debug mode ([PR #173](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/173)).
- Fixed string matching for EASE tile file to accommodate new "EASE*-Pfafstetter" tile file for runoff routing purposes ([PR #160](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/160)).
- Fixed GEOSlandpert build when MKL is unavailable by enabling MKL-specific code paths only when MKL is detected ([PR #162](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/162)).
- Fixed NAG Fortran compiler issues ([PR #170](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/170)).
- Fixed missing deallocate and nullify statements ([PR #180](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/180)).
  
### Removed

- Removed 2d lfs collection from HISTORY.rc template ([PR #156](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/156)).

-----------------------------

## [v3.2.0] - 2025-11-26

- 0-diff vs. v3.1.0 (except for lat/lon fields in "1d" nc4 output, which have roundoff differences between files directly generated with MAPL [new default] and files generated with tile_bin2nc4 [discontinued]).

### Added

- Added reader for surface meteorological forcing from S2S-3 ([PR #138](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/138)).
- Added matlab reader for binary mwRTM vegopacity file ([PR #142](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/142)).

### Changed

- Changed default format of tile-space HISTORY output to nc4 ([PR #144](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/144)).
- Enable remapping of landice restarts from ldas_setup ([PR #146](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/146)).
- Commented out static QC mask in CYGNSS obs reader ([PR #151](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/151)).
- Cleaned up ldas_setup; split out ldas.py and setup_utils.py; restored ntasks-per-node option ([PR #107](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/107)).
- Update `GEOSlandassim_GridComp/io_hdf5.F90` to allow for use with HDF5 1.14 ([PR #139](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/139)).

### Fixed

- Fixed bug in ASCAT EUMET soil moisture obs reader; bumped max_obs limit ([PR #148](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/148), [PR #151](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/151)).
- Provide default "zoom" value for remap_restarts yaml file ([PR #137](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/137)).
- Fixed Restart=1 when the domain is not global ([PR #107](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/107)).

-----------------------------

## [v3.1.0] - 2025-06-26

- 0-diff vs. v3.0.0.

### Added

- Added python package for post-processing ObsFcstAna output into data assimilation diagnostics ([PR #87](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/87), [PR #111](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/111)).
- Support for 2d output from EASE tile space and 2d output on EASE grid:
  - Switched EASE grid handling to new MAPL EASE Grid Factory ([PR #115](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/115)).
  - Revised pre-processing of HISTORY template ([PR #118](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/118)).
- Support for tile space of stretched cube-sphere grids ([PR #109](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/109)).

### Changed

- Revised experiment setup for coupled land-atm DAS ([PR #102](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/102)).
- Updated defaults in LDASsa_DEFAULT_inputs_*.nml files ([PR #104](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/104)).
- Added optional SLURM "constraint" ([PR #112](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/112)).
- Specify only "ntasks_model" in SLURM resource request ([PR #106](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/106)).

### Fixed

- UDUNITS error ([PR #101](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/101), [PR #123](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/123)).

-----------------------------

## [v3.0.0] - 2025-05-28

- 0-diff vs. v2.0.0.

### Added

- Added functionality to simulate landice tiles ([PR #18](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/18)).
- Added functionality to read nc4-formatted tile file ([PR #18](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/18)).
- Added model-based QC of (MODIS) snow cover area fraction observations using layer-1 soil temperature ([PR #96](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/96)).
- Added default settings and command line args for coupled land-atm DAS ([PR #94](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/94)).

------------------------------

## [v2.0.0] - 2025-04-15

- 0-diff vs. v1.1.0.

### Added

- New update_type for joint 3d soil moisture and 1d snow analysis (Tb+sfmc+sfds+SCF obs) ([PR #68](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/68)).
- Updated subroutine read_obs_sm_ASCAT_EUMET() to work with both original and revised file name templates ([PR #69](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/69)).
- Added CYGNSS soil moisture reader ([PR #76](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/76)).
- Added M21C surface met forcing ([PR #77](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/77)).
- Added Github Actions workflow for testing and building GEOSldas_GridComp ([PR #86](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/86)).

### Changed

- Revised variable names (SHORT_NAME) and descriptions (LONG_NAME) to match M21C file specs ([PR #72](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/72)).
- Updated subroutines read_obs_sm_ASCAT_EUMET(), read_obs_SMAP_halforbit_Tb(), read_obs_SMOS() and read_obs_MODIS_SCF() with hardcoded time ranges for when observations are available and should be read ([PR #73](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/73)).
- Renamed tilecoord%pfaf to %pfaf_index; added matlab tile file reader ([PR #78](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/78)).
- Improved setup of coupled land/atm DAS (incl. changed nomenclature of met forcing files: "Nx+-" --> "bkg.lfo_*") ([PR #81](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/81)).

### Removed

- Removed support for SLES12 operating system at NCCS ([PR #83](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/83)).

-----------------------------

## [v1.1.0] - 2024-11-05

- 0-diff vs. v1.0.2 except for data assimilation in cube-sphere tile space ([PR #41](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/41)).

### Changed

- More optimal distribution of tiles on processors for cubed-sphere tile space ([PR #41](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/41)).
- Updates to scripting to allow for Intel MPI ([PR #57](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/57)).

### Fixed

- Do not increment CO2_YEAR when it is a no-data-value; for Catchment simulations, exclude CatchCN-specific resource variables from LDAS.rc ([PR #51](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/51)).
- Bug fix and improved efficiency in matlab script for generation of mwRTM_param ([PR #46](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/46)).
- Changed EXPDIR to absolute path for POSTPROC_HIST>0 option to work ([PR #42](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/42)).
- Support HISTORY output of ASNOW alone from ENSAVG Gridcomp  ([PR #49](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/49)).

### Removed

- Remove restart options F and G  ([PR #40](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/40)).

-----------------------------

## [v1.0.2] - 2024-04-12

- 0-diff vs. v1.0.1.

### Fixed

- Bug fix for state increment array referencing in update_type=13 ([PR #26](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/26)).
- Fixed CI for LDAS workflow ([PR #34](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/34)).

### Changed

- Moved external `GEOSgcm_GridComp` repository to under `GEOSldas/src/Components` for
  consistency with directory structure of GEOSgcm and GEOSadas  ([PR #27](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/27), [PR #30](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/30)).
- Changed lenkf.j.template to python string ([PR #16](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/16)).


-----------------------------

## [v1.0.1] - 2024-04-10

- 0-diff vs. v1.0.0.

### Fixed

- ldas_setup: Changed entry 'slurm' to 'slurm_pbs' to match remap_params.tpl ([PR #17](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/17)).

-----------------------------

## [v1.0.0] - 2024-03-26

- Inaugural version.  0-diff vs. GEOSldas v18.0.0.

-----------------------------
