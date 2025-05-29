# Changelog for `GEOSldas_GridComp`

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

-----------------------------

## [Unreleased]

### Added

### Changed

- linked EASE_conv from MAPL
- Regridding from EASE to other grid types

### Fixed

### Removed

### Deprecated

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


