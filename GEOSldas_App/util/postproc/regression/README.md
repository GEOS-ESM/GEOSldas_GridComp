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

Run the standard start/stop tests:

    cd /discover/nobackup/.../EXPID
    ./regress/start_stop_model.sh

To include the optional layout-invariance test (T4):

    RUN_LAYOUT=1 ALT_1D=120 ./regress/start_stop_model.sh

See the “Layout-Invariance Test (T4)” section below for details.


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
    
# Layout-Invariance Test (T4)

In addition to the core start/stop regression (T1–T3), the script supports an
optional **layout-invariance test (T4)**.  
This test verifies that GEOSldas produces identical results when the number of
MPI tasks along the active axis (NX or NY) is changed.

T4 is disabled by default.  

Enable it by running:

    RUN_LAYOUT=1 ALT_1D=<npes> ./regress/start_stop_model.sh

where `ALT_1D` is the alternate number of MPI tasks (e.g., 84, 120, 126, ...).

### What T4 does

- Uses the T1 run directory as a frozen template.
- Creates a new sandbox sub-experiment (`run_T4/`) with:
  - identical model configuration,
  - different number of MPI tasks (`ALT_1D`),
  - identical tile distribution file (IMS/JMS.rc) pre-built by `preprocess_ldas.x`.

- Runs a 24-hour GEOSldas simulation under the alternate layout.
- Compares **T1 (baseline)** vs **T4 (alternate layout)**:

  **HISTORY (tolerant compare):**
  - Uses `nccmp -dmfMNS -G history -t <ABS_TOL> -T <REL_TOL>`  
    allowing tiny floating-point differences from MPI reduction ordering.

  **Restarts (strict compare):**
  - Uses `nccmp -dmfgMNS`  
    and requires bit-for-bit identical restart fields at the final time.

### When T4 passes

A passing T4 test means:

- Changing layout (task decomposition) does **not** affect model results,
- HISTORY fields agree within tolerance,
- Restart fields agree bit-for-bit.

### When T4 is skipped

T4 runs only if `RUN_LAYOUT=1` is provided.  
Normal users running only T1–T3 do not trigger layout testing.


# Environment variables

| Variable              | Description                               | Default                                 |
|-----------------------|-------------------------------------------|-----------------------------------------|
| EXPDIR                | Experiment root (run/, input/, build/, output/) | auto-detected                     |
| EXPDOMAIN             | Domain under output/                       | auto-detected                          |
| SUBMIT                | Batch command (Slurm only)                 | `sbatch`                               |
| ALT_1D                | Alternate 1-D task count for layout test   | required if `RUN_LAYOUT=1`             |
| ABS_TOL               | Absolute tolerance for tolerant nccmp      | `1e-15`                                |
| REL_TOL               | Relative tolerance for tolerant nccmp      | `1e-12`                                |
| NCCMP_FLAGS_TOL       | Tolerant compare flags                     | `-dmfMNS -G history -t 1e-15 -T 1e-12` |
| HIST_STEP_SEC         | Step for HISTORY collection                | 21600 (6 h)                            |
| HIST_STEP_OFFSET_SEC  | Center offset (+3 h)                       | 10800                                  |


# Comparison logic

GEOSldas start/stop regression uses three comparison modes:

### 1. Restart files (strict compare)
Restart files are compared with full data + metadata + attributes:

    nccmp -dmfgMNS fileA fileB

This must be bit-for-bit identical for the test to pass.

### 2. HISTORY files (data-only strict compare)
For HISTORY collections, only variables are compared (metadata ignored):

    nccmp -dNM fileA fileB

### 3. Layout-invariance tests (tolerant compare)
For layout (T4) or tolerant mode, the script uses:

    nccmp -dmfMNS -G history -t <ABS_TOL> -T <REL_TOL>

which is controlled by:

    NCCMP_FLAGS_TOL = -dmfMNS -G history -t 1e-15 -T 1e-12

This tolerates tiny floating-point differences caused by MPI layout changes.


# Notes

- The 6-hour profile is used for both CF (2d/Nx) and EASE (1d/Nt).  
  It reduces runtime and I/O while staying bit-for-bit safe for segmented runs.


# Maintenance

Templates (`templates/HISTORY_1d.rc`, `templates/HISTORY_2d.rc`) are version-controlled.  
If land-ice is disabled, the `glc` stream is ignored automatically by GEOSldas.
