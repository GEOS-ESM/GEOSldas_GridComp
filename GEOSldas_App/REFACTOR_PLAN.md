# ldas_setup refactor: plan and status

Scope: `ldas_setup`, `ldas.py`, `setup_utils.py`, `lenkf_j_template.py` (GEOSldas_App).
This file is the single place to track progress. Update the **Status** table and the
**Log** at the end of each work session.

Guiding rule: the experiment that gets created (rc files, `lenkf.j`, `ldas_batchrun.j`,
symlinks, copied restarts/BCs) must stay **identical** to the reference unless the
change is an intentional bug fix listed below.

**Regression check (run after every change):**
```
cd GEOSldas_App
python3 -m unittest discover -s tests -v          # unit tests
python3 tests/mock_regression.py                  # end-to-end, ref commit vs working tree
```
`mock_regression.py` compares against `REF_COMMIT` (0f6155f, pre-refactor) by default.
For an intentional behavior change (bug fix), run it, check that the diff matches what you
expect, and record that diff in the bug entry below.

---------------------------------------------------------------------------------------

## Status

| Phase | Description                                    | Status        | Commit |
|-------|------------------------------------------------|---------------|--------|
| 1     | Add helpers to setup_utils                     | DONE          | 0f6155f |
| 0     | Safety net: unit tests, diff tool, mock e2e harness | DONE (Discover run still recommended) | see Log |
| 1b    | Use the helpers in ldas.py / ldas_setup        | DONE          | see Log |
| 2     | Fix the known bugs (B1-B3, B5, B6, B8-B10; B4/B7 are not bugs) | DONE | see Log |
| 3     | Split `ldas.__init__`                          | TODO          | |
| 4     | Table-driven restart staging (`createLnRstBc`) | TODO          | |
| 5     | rc file generation (`createRCFiles`)           | TODO          | |
| 6     | Job scripts, site config, small helpers        | TODO          | |
| 7     | (optional) Split into modules                  | TODO          | |

### Phase 0 details (DONE)
- `tests/test_setup_utils.py`: stdlib `unittest` (pytest is not available locally; python 3.6).
  Covers parseInputFile (comments, `GEOSldas=>`/`GEOSagcm=>` prefixes, duplicate keys, bad lines,
  missing file), require, glob_first, bcs_version_from_path, hours_to_hhmmss, replace_in_file,
  read_gridname (.til and .til.domain), and the key lists.
- `tests/compare_exp.py REF_EXPDIR NEW_EXPDIR [--same P_REF P_NEW]...`: diffs two experiment
  trees. It compares the layout, symlink targets, text in run/ and BCS_info.txt (with exp paths
  and `--same` prefixes normalized), and binary files byte-for-byte.
- `tests/mock_regression.py`: end-to-end harness that needs no Discover and no build.
  - It builds a fake install from `git show REF_COMMIT:` and another from the working tree.
    The `@GEOS_SITE@` placeholder is set to NCCS. The stubs are preprocess_ldas.x (python),
    GEOSldas.x, and the remap_* modules. The real pieces are ldas_setup, ldas.py,
    setup_utils.py, lenkf_j_template.py, process_hist.csh (csh), GEOSldas_*.rc, nml files,
    and GEOS_SurfaceGridComp.rc (from ../@GEOSgcm_GridComp).
  - It creates fake BCS, restart experiment, MWRTM and met data, runs `ldas_setup setup` with
    both installs, and compares the results.
  - 14 cases: r1_global, r1_zoom, r1_ens_assim (nens=2, LAND_ASSIM, PERTURBATIONS, MWRTM,
    ENSEMBLE_FORCING), r1_zoom_ens_assim, r3_catchcn, r1_lake_landice_route,
    r1_lake_landice_zoom, r1_issm, r2_remap, rM_merra2, multi_sgmt, job_sgmt_days,
    nml_hist_custom (NML_INPUT_PATH, POSTPROC_HIST, irrig, aerosol, all optional batinp keys),
    tile_txt.
  - Options: `--case NAME..`, `--keep` (keep the work dir to inspect it), `-v`, `--ref COMMIT`.
  - Added in Phase 2: coupled mode `ladas_cpl1`/`ladas_cpl2` (the /discover paths are replaced
    in the fake installs), the NAS site, and one case per bug fix. Cases have an
    `expect` field: same | diff | ref_fails | new_fails.
  - NOT covered: real Fortran zoom/optimize and real remapping.
- Recommended but not done: run a few real cases on Discover (incl. `--ladas_cpl 1/2`) with
  installs of the old and new code, and compare them with `tests/compare_exp.py`.

### Phase 1b details (DONE)
- `ldas.py`:
  - every `assert` changed to `require(...)` (asserts are skipped under `python -O` and print tracebacks)
  - every `exit("...")` changed to `raise LdasSetupError(...)`
  - all `print cmd` + `sp.call(shlex.split(cmd))` pairs changed to `run_cmd(cmd)` (10 places)
  - `glob(...)[0]` changed to `glob_first(pattern, what)` for the BC files, vegdyn, and mk_restarts restarts
  - added `require` checks before indexing the tile-file lists (RESTART=2, TILE_FILE_FORMAT=TXT)
  - the 3 BCS_VERSION loops changed to `bcs_version_from_path`
  - `fileinput` in-place edits changed to `replace_in_file` (same output: rstrip + replace)
  - removed the unused `fileinput` and `shlex` imports
  - the coupled-mode argument checks are now a loop that names the missing `--arg`
- `ldas_setup`: added `main()`, removed the `status`/`assert(status)` pattern, and it now catches
  `LdasSetupError` -> prints `ldas_setup ERROR: <msg>` to stderr and exits with 1.
- Behavior changes: only the error messages and exit codes. mock_regression: 14/14 IDENTICAL.

---------------------------------------------------------------------------------------

## Known bugs (Phase 2): DONE

Each fix is checked by a dedicated case in `tests/mock_regression.py` (shown in brackets).

- **B1 (FIXED)** A missing DO_ISSM crashed setup: TypeError via `int(ExeInputs.get('DO_ISSM'))` when
  landice tiles were used, and KeyError in the process_hist.csh call in *every* run.
  **This also broke the coupled land-atm DAS (`--ladas_cpl 1/2`) completely**, because the plain
  HIST.rc is always processed and coupled runs have no exeinp.
  Fix: `self.do_issm = int(ExeInputs.get('DO_ISSM', 0))`, which is used for with_issm and
  process_hist. It is not written to LDAS.rc, so LDAS.rc is unchanged.
  [B1_no_do_issm, ladas_cpl1, ladas_cpl2: ref fails, new OK]
- **B2 (FIXED)** MWRTM_PATH without mwRTM_param.nc4: the old check could never pass (so it always
  aborted with a confusing message), and the `del` after it was dead code. Kept the
  behavior (abort) with a clear message: `MWRTM_PATH/BCS_RESOLUTION [..] must contain mwRTM_param.nc4`.
  [B2_mwrtm_missing: new fails]
- **B3 (FIXED)** RUN_ROUTE with RESTART=M, and DO_ISSM=1 with RESTART=M, are now rejected in
  `__init__` before any remapping. The later (now unreachable) checks were removed. [B3_route_merra2]
- ~~**B4**~~ **NOT A BUG (intended, confirmed 2026-10-22).** If the exeinp file exists,
  `ladas_cpl` is forced to 0 and `--ladas_cpl` is ignored. Keep this behavior in all phases.
  Done: removed the duplicated `int(self.ladas_cpl)`.
- **B5 (FIXED)** LSM_CHOICE now goes through the `LSM_CHOICE_TO_CATCH` table. An unknown value
  -> error (it used to leave `catch=''`). The existing "3 not supported with land tiles"
  check is unchanged. [B5_bad_lsm_choice]
- **B6 (FIXED)** `ldas_batchrun.j` on NAS now uses `jobidN=$(qsub [-W depend=afterok:$jobidN-1] lenkf.j)`
  (it used to use sbatch syntax). NCCS output is unchanged. [B6_nas_multi_sgmt, nas_single: diff in
  ldas_batchrun.j only. NOTE: the qsub syntax has not been tested on a real NAS system]
- ~~**B7**~~ **NOT A BUG (intended, confirmed 2026-10-22).** `_calculateJobSegments` uses
  month = 30 days and year = 365 days on purpose. Do NOT switch to `relativedelta`, and keep the
  segment dates byte-identical in all phases. Done: removed the debug `print`s (console only).
  (Side effect: with the default JOB_SGMT a run from Jan 1 to Feb 1 gets 2 segments.)
- **B8 (FIXED, soft)** `run_cmd(cmd, check=False)` now prints a WARNING when the return code is
  non-zero, and raises if `check=True`. All calls still use the default (warn only), because the
  exit codes of preprocess_ldas.x / process_hist.csh on Discover have not been confirmed.
  TODO (later, on Discover): switch preprocess_ldas.x calls to `check=True`.
- **B9 (FIXED)** The ISSM_INTERNAL_CHECKPOINT_FILE key is now set only on the first loop iteration,
  which keeps the LDAS.rc key order. (The value is the same as the generic per-key one; the order
  is fixed properly in Phase 5.)
- **B10 (FIXED)** Constants at the top of ldas.py: `LADAS_BCS_PATH`, `LADAS_MWRTM_PATH`, `ROUTE_RST_DIR`
  (the mock harness redirects the first two to fake data for the coupled cases).

Other cleanups to do along the way: temp files (`domain_def`, f2g, optimize output) are never
deleted; `os.chdir(self.rundir)` side effects in createRunScripts/createBatchRun.

---------------------------------------------------------------------------------------

## Phases 3-7 (design)

**Phase 3: split `__init__` (~500 lines)** into ordered private steps:
`_read_cmdline()`, `_apply_ladas_coupling()`, `_set_defaults()`,
`_set_model_options()` (LSM/catch, tile types, ensemble ids), `_calculateJobSegments()`,
`_resolve_bcs_files()` (tile + BC globbing), `_check_restart_inputs()`,
`_read_resource_inputs()`, `_set_paths()`. Move the domain-namelist temp file into
`createLnRstBc`. Put the coupled-mode constants in one block at the top of the module.

**Phase 4: `createLnRstBc` (~490 lines)**
- `_link_met_forcing()`, `_stage_tile_file()`, `_compute_zoom()` (f2g -> zoom-by-type dict),
  `_stage_bcs()`, `_link_issm_bcs()`, `_copy_landassim_seeds()`, `_remap_restarts()`,
  `_stage_restarts()`, `_stage_mwrtm()`, `_write_bcs_info()`.
- Component table that drives `_stage_restarts`, roughly:
  `name, enabled, all_members?, src(RESTART 1/3), mk_restarts pattern(2/M), zoom subcmd, link name`
  for catch/vegdyn/lake/landice/issm/route/landpert. One generic
  `locate -> zoom-or-copy -> link` routine. Keep these special cases:
  members after 0 fall back to member-0 files (catch/vegdyn/route), lake/landice/issm only
  for iens==0, and vegdyn is not under Yyyyy/Mmm.
- This phase is the highest risk: run **all** baseline cases.

**Phase 5: `createRCFiles`**
- `_archive_inputs()` (copy-with-suffix helper), `_copy_nml()`, `_optimize_layout()`,
  `_write_history_rc()`, `_write_cap_rc()`, `_write_ldas_rc()`.
- `_build_ldas_rc_dict()` should be a pure function (no file I/O) with unit tests.

**Phase 6: job scripts and helpers**
- Date helpers `y4m2(dt)` -> 'Yyyyy/Mmm', `stamp(dt)` -> 'yyyymmdd_hhmm' (used 5+ times).
- `SITE = {'NCCS': {...qos, constraint, submit, dependency syntax}, 'NAS': {...}}`.
- `createBatchRun` uses the site config (fixes B6 if not already done).

**Phase 7 (optional): modules**: `ldas_config.py`, `ldas_restarts.py`, `ldas_rc.py`,
with `ldas.py` as the orchestrator. Remember to:
- add new .py files to `scripts` in `CMakeLists.txt`
- only `ldas.py` goes through `configure_file` (`@GEOS_SITE@`, `@MPI_STACK@`). Keep those
  placeholders in a single small module and pass them to the other modules.

---------------------------------------------------------------------------------------

## How to resume

1. `git status` / `git log --oneline -5` in `@GEOSldas_GridComp`. Compare with the Log below.
2. `cd GEOSldas_App && python3 -m unittest discover -s tests && python3 tests/mock_regression.py`
   (both should pass).
3. Pick the next TODO in the Status table. Use one commit per bug or sub-step.
4. Optional, on Discover: build/install and compare real cases with `tests/compare_exp.py`.
5. Update this file (Status, phase details, Log).

### Local smoke test (no Discover; mock_regression.py is the more complete option)
```
S=/tmp/stubs; mkdir -p $S
printf 'STRETCH_GRID={}\ndef init_merra2(a): pass\ndef get_config_from_answers(*a,**k): return {}\n' > $S/remap_utils.py
echo 'other_restarts=None' > $S/remap_other_restarts.py
echo 'catchANDcn=None'     > $S/remap_catchANDcn.py
PYTHONPATH=$S python3 ldas_setup sample --batinp
PYTHONPATH=$S python3 ldas_setup setup /nonexistent x y      # -> "ldas_setup ERROR: ...", rc=1
```

---------------------------------------------------------------------------------------

## Log

- 2026-10-20: Phase 1 committed (0f6155f): helpers in setup_utils, GRIDNAME fix.
- 2026-10-22: Review + plan written. Phase 0 tooling (unit tests, compare_exp.py) and
  Phase 1b done in the working tree (not committed). Next: run baselines on Discover and
  compare. If there are 0 diffs, commit Phase 0 + 1b, then start Phase 2 with B1.
- 2026-10-22: User confirmed that B4 (exeinp present -> ladas_cpl=0) and B7 (30-day month /
  365-day year in job segments) are intended behavior. Removed them from the Phase 2 fix list.
  IDs are not renumbered.
- 2026-10-22: Added tests/mock_regression.py (end-to-end mock harness, 14 cases). Phase 1b
  passes 14/14 IDENTICAL against 0f6155f. Committed Phase 0 + 1b.
- 2026-10-22: Phase 2 done (B1-B3, B5, B6, B8-B10 + B4/B7 cleanups). Harness: 22 cases, PASS.
  Important: B1 means coupled mode crashed in 0f6155f.
