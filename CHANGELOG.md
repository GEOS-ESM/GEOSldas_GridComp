# Changelog for `GEOSldas_GridComp`

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

-----------------------------

## [Unreleased]

### Added

- Added landice grid comp
- New update_type for joint 3d soil moisture and 1d snow analysis (Tb+sfmc+sfds+SCF obs)

### Changed

- Updated subroutine read_obs_sm_ASCAT_EUMET() to work with both original and revised file name templates. 
- Updated subroutines read_obs_sm_ASCAT_EUMET(), read_obs_SMAP_halforbit_Tb(), read_obs_SMOS() and read_obs_MODIS_SCF() with hardcoded time ranges for when observations are available and should be read.
- Revised variable names (SHORT_NAME) and descriptions (LONG_NAME) to match M21C file specs.

### Fixed

### Removed

### Deprecated

-----------------------------

## [v1.1.0] - 2024-11-05

- 0-diff vs. v1.0.2 except for data assimilation in cube-sphere tile space ([PR #41](https://github.com/GEOS-ESM/GEOSldas_GridComp/pull/41)).

### Changed

- More optimal distribution of tiles on processors for cubed-sphere tile space.
- Updates to scripting to allow for Intel MPI.

### Fixed

- Do not increment CO2_YEAR when it is a no-data-value. For Catchment simulations, exclude CatchCN-specific resource variables from LDAS.rc.
- Bug fix and improved efficiency in matlab script for generation of mwRTM_param.
- Changed EXPDIR to absolute path for POSTPROC_HIST>0 option to work.
- Support HISTORY output of ASNOW alone from ENSAVG Gridcomp.

### Removed

- Remove restart options F and G.

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


