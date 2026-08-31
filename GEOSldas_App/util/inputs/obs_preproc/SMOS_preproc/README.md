# SMOS_preproc

## 1. Purpose
Preprocesses SMOS L1C brightness temperature data for near-real-time (NRT)
ingestion into the MERRA21C-land pipeline:

```
EE .zip  --(smos-ee-to-nc.sh)-->  NetCDF  --(preprocess_nc)-->  EASEv2 M36 "REG" binaries  --(SCLF1C_reg2fit)-->  Tb40 "FIT" binaries
```

Normal operation is a **daily cron job** that processes the previous day's
data. It also supports on-demand backfill of a single day or a date range.

## 2. Requirements
- Runs on Discover/NCCS.
- Environment: `module load python/GEOSpyD` (provides numpy, scipy, netCDF4,
  PyYAML — no other Python environment should be needed).
- External script dependency: `smos-ee-to-nc.sh` (path set in `config.yaml`,
  currently maintained under `/discover/nobackup/projects/gmao/smap/...` —
  **not part of this repo**, 
- Static input: `data/GEOSIT_to_EASEv2_M36.mat` — pre-generated EASEv2
  regridding weights, link to the actual location on Discover
- Static Aux files: `data/SM_OPER_AUX_GAL_SM_20050101T000000_20500101T000000_001_003_3`
  link to the actual location on Discover. 

## 3. Configuration (`config.yaml`)
| Key                         | Meaning                                                                                         |
|-----------------------------|-------------------------------------------------------------------------------------------------|
| `ee_to_nc_script`           | Path to the external `smos-ee-to-nc.sh` converter.                                              |
| `smos_base_path`            | Where incoming SMOS `SM_*_MIR_SCLF1C_*.zip` (EE) files land, organized `<base>/Y<yyyy>/M<mm>/`. |
| `tmp_nc_path`               | Scratch directory for converted NetCDF files. A `to_delete/` subfolder is created here. Files therein are safe to delete after completion. |
| `out_reg_path`              | Output root for REG binaries (`SMOS_reg_Tb_*.bin`, organized by `<out_reg_path>/<YYYYMM>/`). FIT binaries are written to a sibling directory: `_reg_` in this path is replaced by `_fit_` and further nested under `SMOS_fit_poly2/<YYYYMM>/`. |
| `GEOSIT_path`               | Path to GEOS-IT data.                                                                           |
| `GEOSIT_to_EASEv2_M36_file` | Path and name of file mapping from GEOS-IT output grid to EASEv2_M36 grid.                      |
| `SM_OPER_AUS_GAL_SM_path`   | Path to SMOS galaxy correction files.                                                           |


Current values are set for the production Discover paths under
`/discover/nobackup/projects/gmao/smap/SMAP_Nature/SMOS/...` and
`/discover/nobackup/dao_ops/SMOS/AOSMOS.4662/SCLF1C/`.

## 4. Usage
```
module load python/GEOSpyD

python SMOSproc_main.py                                    # NRT: process yesterday's data (normal cron mode)
python SMOSproc_main.py --date  20260801                   # Backfill a single day
python SMOSproc_main.py --start 20260801 --end 20260803    # Backfill an inclusive date range
```
`--date` and `--start`/`--end` are mutually exclusive; `--end` requires
`--start`.

## 5. Layout
```
SMOSproc_main.py          top-level driver 
config.yaml               paths 
data/                     static regridding weights (GEOSIT_to_EASEv2_M36.mat) and AUX 
src/preprocess_nc.py      NetCDF -> EASEv2 M36 REG binaries
src/SCLF1C_reg2fit.py     REG -> Tb40 FIT binaries
src/readwrite/            I/O helpers (SMOS NetCDF, REG binary, GEOS-IT, aux Gal SM)
src/helper/               grid/geometry/time utilities (EASEv2 indexing, h/v tile
                          conversion, celestial angle calc, pentad/leap-year,
                          galactic+atmospheric correction)
```


