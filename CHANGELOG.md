# Changelog for `GEOSldas_GridComp`

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

-----------------------------

## [Unreleased]

### Changed

- Moved the location mepo clones `GEOSgcm_GridComp` to under `GEOSldas/src/Components` parallel
  to `GEOSldas_GridComp` rather than inside `GEOSldas_GridComp` itself. This is to be
  consistent with the structure of the GEOSgcm and GEOSldas components.

-----------------------------

## [v1.0.2] - 2024-04-11

### Fixed

- Bugfix for state increment array referencing in update_type=13

-----------------------------

## [v1.0.1] - 2024-04-10

### Fixed

- ldas_setup: Changed entry 'slurm' to 'slurm_pbs' to match remap_params.tpl

-----------------------------

## [v1.0.0] - 2024-03-26

- Inaugural version.  0-diff vs. GEOSldas v18.0.0.

-----------------------------

### Added

### Changed

### Fixed

### Removed

### Deprecated

