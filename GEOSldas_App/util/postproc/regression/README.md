# GEOSldas Global Regression: Model Start/Stop & Layout

## Overview

This regression is run after you have already built and executed a working GEOSldas experiment.

You must have:

- A complete experiment directory containing:
  - run/
  - input/
  - build/
  - output/<DOMAIN>/
- Valid restart files under:
  output/<DOMAIN>/rs/ens0000/
  (e.g., CURRENT.catch_internal_rst.*, CURRENT.landice_internal_rst.*)
- A run/LDAS.rc that defines your grid type (CF or EASE)

The regression does not modify your experiment.
It makes a self-contained sandbox copy, runs start/stop tests, and compares results.

## What the Regression Does

The regression runs GEOSldas in an isolated sandbox cloned from your experiment.
It forces a 6-hour time-averaged HISTORY profile (small and fast) and verifies that:

1. Restart files are identical between:
   - a 24-hour continuous run, and
   - a 12 h + 12 h split run.

2. HISTORY output (6-hour centers) is identical for the same 24-hour window.

## Grid Handling

This regression is grid-agnostic.

### CF (Cubed Sphere 2-D Grids)
- HISTORY collection: tavg24_2d_*_Nx

### EASE (1-D Grids)
- HISTORY collection: tavg24_1d_*_Nt

Both output types are normalized to:
- 6-hour frequency: 060000
- Reference time: 000000

## Safety

Your real experiment is not modified.

All regression work occurs in:

    regress/sandbox/<EXPID>  


## Regression package layout

<pre>
util/postproc/regression/
├─ start_stop_model.sh        # regression driver
├─ templates/
│   ├─ HISTORY_2d.rc          # CF (2d/Nx) 6-hour tavg only
│   └─ HISTORY_1d.rc          # EASE (1d/Nt) 6-hour tavg only
├─ README.md                  # this file

When a regression run starts, this structure appears under your experiment:

&lt;EXPID&gt;/
├─ run/                       # original job files (unchanged)
├─ input/                     # restart, tile, forcing, etc.
├─ build/                     # model binaries
├─ output/&lt;DOMAIN&gt;/           # real experiment outputs
│   ├─ rs/ens0000/            # restarts (catch, land-ice)
│   ├─ cat/ens0000/           # HISTORY (tavg24_*.nc4)
│   └─ rc_out/                # category files
└─ regress/
    ├─ logs/                  # regression stdout/stderr with timestamps
    ├─ sets/                  # collected results per segment:
    │   ├─ T1_*               # 24 h run
    │   ├─ T2_*               # 12 h first half
    │   └─ T3_*               # 12 h second half
    └─ sandbox/&lt;EXPID&gt;/       # isolated copy used for the run
        ├─ run/               # patched job/rc files
        ├─ build/             # symlink to ../build
        ├─ output/&lt;DOMAIN&gt;/   # new outputs written here
        └─ scratch/           # Slurm log/stdout/err for sandbox runs
</pre>


To inspect the sandbox after a run, comment out the final cleanup line
in `start_stop_model.sh`. By default, the sandbox is deleted after a PASS.


# Quick start

Run your experiment once so that restart files and outputs exist.
The regression uses these restarts as inputs.

Run only start/stop tests in the regression driver:

    cd /discover/nobackup/.../EXPID
    ./regress/start_stop_model.sh

Run with layout test:

    ./regress/RUN_LAYOUT=1 ALT_1D=120 ./start_stop_model.sh

where `ALT_1D` can be 84, 120, 126, etc., depending on grid resolution.


# What the regression does

- Creates `regress/sandbox/<EXPID>` and copies your run directory.
- Detects grid type (CF or EASE) and applies the correct 6-hour HISTORY template.
- Adjusts environment variables:

      DO_HISTORY=TRUE
      DO_HIST=TRUE
      POSTPROC_HIST=0

- Runs:
  - **T1** – single 24-hour job  
  - **T2** – 12-hour run to mid-time  
  - **T3** – 12-hour run to final time

- Compares:
  - **RESTARTS:** T1 (24 h) vs T3 (12 h + 12 h)
  - **HISTORY:** T1 vs [T2 ∪ T3] at 03/09/15/21Z centers


# Environment variables

| Variable              | Description                               | Default                      |
|----------------------|-------------------------------------------|------------------------------|
| EXPDIR               | Experiment root (run/, input/, build/, output/) | auto-detected          |
| EXPDOMAIN            | Domain under output/                       | auto-detected                |
| SUBMIT               | Batch command (Slurm only)                 | `sbatch`                     |
| ALT_1D               | Alternate 1-D task count for layout test   | required if `RUN_LAYOUT=1`   |
| NCCMP_FLAGS_TOL      | Tolerant compare flags                     | `-dmfgqMNS -t 1e-12 -T 1e-6` |
| HIST_STEP_SEC        | Step for HISTORY collection                | 21600 (6 h)                  |
| HIST_STEP_OFFSET_SEC | Center offset (+3 h)                       | 10800                        |


# Comparison logic

- Restarts are compared with:

      nccmp -dmfgqMNS

  If strict compare fails, the script performs a tolerant comparison.

- HISTORY compares all 6-hour stamps in the same 24-hour window.


# Notes

- The 6-hour profile is used for both CF (2d/Nx) and EASE (1d/Nt).  
  It reduces runtime and I/O while staying bit-for-bit safe for segmented runs.


# Maintenance

Templates (`templates/HISTORY_1d.rc`, `templates/HISTORY_2d.rc`) are version-controlled.  
If land-ice is disabled, the `glc` stream is ignored automatically by GEOSldas.
