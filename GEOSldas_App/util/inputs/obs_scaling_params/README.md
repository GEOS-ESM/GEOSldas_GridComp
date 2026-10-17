# Observation scaling parameters

Scripts that compute the climatological mean and standard deviation of
observations and model forecasts from the `ldas_ObsFcstAna` files of a GEOSldas
experiment (usually an open loop).  GEOSldas uses these statistics for the
z-score scaling of observations before assimilation:

```text
obs_scaled = m_mean + (obs - o_mean) * m_std / o_std
```

There is a MATLAB and a Python version.  Both read nc4 ObsFcstAna files and give
the same results (checked on a CF0360 open loop, Oct 2026: lat/lon H SAF ASCAT
and CYGNSS identical up to rounding <1e-14; Tb SMAP/SMOS identical except for
one float32 ulp in a few m_std values).

## Python

| File | Contents |
|---|---|
| `run_obs_scaling_params.py` | command-line driver (`--method latlon` or `--method tb`) |
| `get_obs_scaling_latlon.py` | soil moisture obs (e.g. ASCAT, CYGNSS) on a global lat/lon grid; nc4 output for `scale_obs_sfmc_zscore` |
| `get_obs_scaling_tb.py` | Tb obs (e.g. SMAP, SMOS) in tile space, one orbit per run; binary output for `scale_obs_Tb_zscore` |
| `obs_scaling_utils.py` | helpers for the moving-window loop |

The GEOSldas output files are read with `../../shared/python/read_GEOSldas.py`,
and the EASEv2 transforms are from `../../shared/python/EASEv2.py`.  The scripts
can be run from any directory (`python run_obs_scaling_params.py --help`).

For each day of the climatological year, the statistics are computed from all
obs/forecast pairs (both valid) within a moving window of `--window-days` days
(default 75) over all years of the period, and stored at the center of the
window.  `--window-days` must be an odd number of pentads.  Cells/tiles with
fewer than `--ndata-min` obs (default 20) get no-data.  `--run-months` (default
1-12 plus the months needed to wrap the window) must be consecutive months; 12
may be followed by 1.
Output goes to `<exp-path>/<exp-run>/output/<domain>/stats/<out-dir>/`
(`--out-dir` default `z_score_clim`).  The `ldas_obsparam` file is read from
`rc_out/Y<start>/M<start>/` unless `--obsparam` is given.  The run stops with an
error if there are no ObsFcstAna files for the period.

### `--method latlon`

The grid cell of each obs is computed from its lat/lon with `ceil()` in single
precision, as in GEOSldas (`scale_obs_sfmc_zscore`).  By default all species
given with `--species` are combined (`_sp_ALL_` in the file name); use
`--no-combine-species` for one file per species.  Output: one nc4 file with all
73 pentads (`--no-all-pentads` to skip), and optionally one file per pentad
(`--each-pentad`) or per DOY (`--each-doy`).  `--exclude-tiles-csv` drops obs at
tiles flagged in a CSV file (columns `tilenum` and, by default, `reject`).
`--grid-resolution` (default 0.25 deg) must divide 180 and 360.

Memory: the moving-window sums and the all-pentads output are held in memory for
the whole global grid.  At 0.25 deg with the default 75-day window this is about
8 GiB for the combined species, and about 8 GiB per species with
`--no-combine-species`.

```bash
python run_obs_scaling_params.py --method latlon \
  --exp-path /path/to/exps --exp-run OL1 --domain CF0360x6C_GLOBAL \
  --start 2019-02 --end 2020-08 \
  --species ASCAT_HSAF_META_SM,ASCAT_HSAF_METB_SM,ASCAT_HSAF_METC_SM \
  --prefix CF0360_OL1_HSAF_zscore_stats_
```

### `--method tb`

Selects one H-pol and one V-pol Tb species per angle (`--angles`, default 40)
for the orbit (`--orbit A|D`) whose descr contains `--descr-contains` (default
`SMAP_L1C`).  Statistics are computed on EASEv2 M36 "administering" tiles, or on
all model tiles with `--no-convert-grid`.  As in the MATLAB version, a 95-day
window also requires at least `--ndata-min`/2.5 obs in the 31 days at its
center.  One file per pentad (and optionally per DOY) in the big-endian Fortran
sequential format that GEOSldas and `read_Tb_scale_param.py` read.

```bash
python run_obs_scaling_params.py --method tb \
  --exp-path /path/to/exps --exp-run OL1 --domain SMAP_EASEv2_M36_GLOBAL \
  --start 2015-04 --end 2021-03 --descr-contains SMAP_L1C --orbit D \
  --prefix OL1_SMAP_zscore_stats_
```

## MATLAB

- `Run_get_model_and_obs_clim_stats_latlon_grid.m` with
  `get_model_and_obs_clim_stats_latlon_grid.m`: lat/lon method.
- `Run_get_L4_Tb_scale_SMAP.m` with `get_model_and_obs_clim_stats.m`: Tb method.

The MATLAB scripts use the readers in `../../shared/matlab/`.  Settings are
edited in the `Run_*.m` scripts.
