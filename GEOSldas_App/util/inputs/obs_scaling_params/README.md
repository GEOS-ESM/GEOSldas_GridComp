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

## Method

Both versions use the same method.

- **Climatological year.**  All years of the period are mapped onto one non-leap
  year (2014): obs from, e.g., 15 March 2019 and 15 March 2020 go into the same
  day.  29 February is never read.  Months that are covered by more years of the
  period contribute more obs.
- **Moving window.**  For each day, the sums of obs, obs^2, fcst, fcst^2 and the
  number of obs are kept for each grid cell (lat/lon) or tile (Tb).  The
  statistics of a window of `w_days` days are computed from these sums (the stdv
  is the population stdv, i.e., divided by N), and the window then moves on by
  one day.
- **Obs/fcst pairs.**  Only obs with a valid obs value and a valid fcst value are
  used.  Cells/tiles with fewer than `Ndata_min` obs in the window get no-data.
  For the Tb method with a 95-day window, cells/tiles with fewer than
  `Ndata_min/2.5` obs in the 31 days at the center of the window also get
  no-data.
- **Run months.**  The loop runs through `run_months`, i.e., 1-12 followed by
  the first months again (Python default: 1 to ceil(`w_days`/30); in MATLAB,
  set in the `Run_*.m` scripts).  A window is complete only after `w_days`
  days, so the first pass through the first months only fills the window.  The
  second pass gives the windows that are centered in late December to February.
  `run_months` must be consecutive (12 may be followed by 1).
- **Years.**  Each month is read only for the years in which it falls within the
  period.  For example, with a period from 2019-02 to 2020-08, January is read
  for 2020 only, September to December for 2019 only, and February to August for
  both years.  (Python and `Run_get_model_and_obs_clim_stats_latlon_grid.m`
  compute these years from the period; `Run_get_L4_Tb_scale_SMAP.m` sets them
  for each entry of `run_months`.)
- **DOY of the statistics.**  The statistics of a window are stored under
  DOY(end - `w_days`/2), where "end" is the time of the last ObsFcstAna file of
  the last day of the window.  With 3-hourly files (last file at 21z) the window
  is centered on that DOY.  With one file per day (at 0z) its center is one day
  later.  The pentad files contain the statistics of the DOY at the center of
  the pentad (DOY 3 for pentad 1, DOY 8 for pentad 2, ...).
- **lat/lon grid cells.**  The grid cell of each obs is computed in the same way
  as in GEOSldas (`scale_obs_sfmc_zscore`), so that an obs on a cell boundary is
  counted in the cell in which GEOSldas looks up its scaling parameters.  All
  obs in a cell are used, including several obs in the same cell and file.
- **Tb tiles.**  By default, the statistics are computed for one EASEv2 M36
  "administering" tile per M36 grid cell.  This is the model tile at the
  (slightly shifted) center of the M36 cell, unless the cell center is outside
  the bounding box of that tile and more than 20 km from its center.

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

The window length is `--window-days` (default 75; must be an odd number of
pentads), the minimum number of obs is `--ndata-min` (default 20), and
`--run-months` overrides the default run months (see "Method").  The files are
read one at a time, so the run time grows with the number of ObsFcstAna files.
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
all model tiles with `--no-convert-grid`.  The run stops with an error if an
obs is not on an output tile, if a tile has more than one obs of a species in a
file, or if the lat/lon of an obs is not close to that of its tile.  One file
per pentad (and optionally per DOY) in the big-endian Fortran
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
