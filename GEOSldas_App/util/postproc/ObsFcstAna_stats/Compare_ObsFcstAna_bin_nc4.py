#!/usr/bin/env python3

"""
Find and compare matching ObsFcstAna .bin and .nc4 files.

Usage examples:
  ./Compare_ObsFcstAna_bin_nc4.py /path/to/ana/ens_avg
  ./Compare_ObsFcstAna_bin_nc4.py /path/to/ana/ens_avg --max-pairs 5
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

THIS_DIR = Path(__file__).resolve().parent
sys.path.append(str((THIS_DIR / '../../shared/python').resolve()))
from read_GEOSldas import read_ObsFcstAna  # noqa: E402


INT_MAP = [
    ('obs_assim', 'assim'),
    ('obs_species', 'species'),
    ('obs_tilenum', 'tilenum'),
]

FLOAT_MAP = [
    ('obs_lon', 'lon'),
    ('obs_lat', 'lat'),
    ('obs_obs', 'obs'),
    ('obs_obsvar', 'obsvar'),
    ('obs_fcst', 'fcst'),
    ('obs_fcstvar', 'fcstvar'),
    ('obs_ana', 'ana'),
    ('obs_anavar', 'anavar'),
]

TIME_KEYS = ['year', 'month', 'day', 'hour', 'minute', 'second', 'dofyr', 'pentad']


def _as_array(x):
    return np.asarray(x)


def _normalize_missing(a, fill_value=-9999.0):
    """Replace fill_value sentinels with NaN for consistent comparison.

    The binary reader (read_ObsFcstAna) converts obs_obsvar/fcst/fcstvar/ana/
    anavar to NaN but leaves obs_obs, obs_lon, and obs_lat as raw -9999.
    The NC4 reader uses the per-variable missing_value attribute (also -9999).
    Calling this function on both sides guarantees a uniform NaN-based
    comparison regardless of which fields each reader normalizes internally.
    """
    out = np.asarray(a, dtype=np.float64).copy()
    out[np.isclose(out, fill_value)] = np.nan
    return out


def read_obsfcstana_nc4(fname):
    out = {}
    date_time = {}

    with Dataset(fname, mode='r') as nc:
        for bkey, nkey in INT_MAP:
            out[bkey] = _as_array(nc.variables[nkey][:]).astype(np.int32)

        for bkey, nkey in FLOAT_MAP:
            v = _as_array(nc.variables[nkey][:]).astype(np.float64)
            fillv = getattr(nc.variables[nkey], 'missing_value', -9999.0)
            out[bkey] = _normalize_missing(v, fill_value=float(fillv))

        for key in TIME_KEYS:
            date_time[key] = int(getattr(nc, key, -9999))

        out['N_obsf_attr'] = int(getattr(nc, 'N_obsf', out['obs_species'].size))

    return date_time, out


def compare_arrays_int(a, b):
    if a.shape != b.shape:
        return {'shape_match': False, 'n_bad': max(a.size, b.size), 'max_abs': np.nan}
    bad = (a != b)
    return {
        'shape_match': True,
        'n_bad': int(np.count_nonzero(bad)),
        'max_abs': int(np.max(np.abs(a.astype(np.int64) - b.astype(np.int64)))) if a.size > 0 else 0,
    }


def compare_arrays_float(a, b, rtol, atol):
    if a.shape != b.shape:
        return {'shape_match': False, 'n_bad': max(a.size, b.size), 'max_abs': np.nan}

    bad = ~np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
    n_bad = int(np.count_nonzero(bad))

    diff = np.abs(a - b)
    max_abs = float(np.nanmax(diff)) if diff.size > 0 else 0.0

    return {'shape_match': True, 'n_bad': n_bad, 'max_abs': max_abs}


def find_pairs(root):
    bin_files = sorted(root.rglob('*.ldas_ObsFcstAna.*z.bin'))
    nc4_files = sorted(root.rglob('*.ldas_ObsFcstAna.*z.nc4'))

    pairs = []
    missing_nc4 = []
    for b in bin_files:
        n = b.with_suffix('.nc4')
        if n.exists():
            pairs.append((b, n))
        else:
            missing_nc4.append(b)

    missing_bin = []
    for n in nc4_files:
        b = n.with_suffix('.bin')
        if not b.exists():
            missing_bin.append(n)

    return pairs, missing_nc4, missing_bin


def compare_one(bin_file, nc4_file, rtol, atol):
    b = read_ObsFcstAna(str(bin_file))
    t_nc4, n = read_obsfcstana_nc4(str(nc4_file))

    issues = []
    n_bad_total = 0

    # Time/header checks
    t_bin = b['date_time']
    for key in TIME_KEYS:
        key_bin = 'min' if key == 'minute' else ('sec' if key == 'second' else key)
        vb = int(t_bin.get(key_bin, -9999))
        vn = int(t_nc4.get(key, -9999))
        if vb != vn:
            issues.append(f'time mismatch: {key} bin={vb} nc4={vn}')

    n_obs_bin = int(_as_array(b['obs_species']).size)
    n_obs_nc4 = int(_as_array(n['obs_species']).size)
    if n_obs_bin != n_obs_nc4:
        issues.append(f'N_obs mismatch: bin={n_obs_bin} nc4={n_obs_nc4}')
    if int(n.get('N_obsf_attr', n_obs_nc4)) != n_obs_nc4:
        issues.append(f'N_obsf attr mismatch: attr={n.get("N_obsf_attr")} actual={n_obs_nc4}')

    # Integer fields
    for bkey, _ in INT_MAP:
        cb = compare_arrays_int(_as_array(b[bkey]), _as_array(n[bkey]))
        if not cb['shape_match'] or cb['n_bad'] > 0:
            issues.append(f'{bkey}: n_bad={cb["n_bad"]}, max_abs={cb["max_abs"]}')
            n_bad_total += cb['n_bad']

    # Float fields
    for bkey, _ in FLOAT_MAP:
        a_bin = _normalize_missing(_as_array(b[bkey]))
        a_nc4 = _normalize_missing(_as_array(n[bkey]))
        cb = compare_arrays_float(a_bin, a_nc4, rtol=rtol, atol=atol)
        if not cb['shape_match'] or cb['n_bad'] > 0:
            issues.append(f'{bkey}: n_bad={cb["n_bad"]}, max_abs={cb["max_abs"]:.6g}')
            n_bad_total += cb['n_bad']

    return issues, n_bad_total, n_obs_bin


def main():
    parser = argparse.ArgumentParser(description='Compare GEOSldas ObsFcstAna .bin and .nc4 files.')
    parser.add_argument('root', nargs='?', default='.', help='Root directory to scan recursively.')
    parser.add_argument('--rtol', type=float, default=0.0, help='Relative tolerance for float compare.')
    parser.add_argument('--atol', type=float, default=1e-6, help='Absolute tolerance for float compare.')
    parser.add_argument('--max-pairs', type=int, default=0, help='Max number of pairs to compare (0=all).')
    args = parser.parse_args()

    root = Path(args.root).expanduser().resolve()
    pairs, missing_nc4, missing_bin = find_pairs(root)

    print(f'scan_root: {root}')
    print(f'found pairs: {len(pairs)}')
    print(f'missing nc4 for bin: {len(missing_nc4)}')
    print(f'missing bin for nc4: {len(missing_bin)}')

    if len(missing_nc4) > 0:
        print('example missing nc4:')
        for p in missing_nc4[:5]:
            print(f'  {p}')

    if len(missing_bin) > 0:
        print('example missing bin:')
        for p in missing_bin[:5]:
            print(f'  {p}')

    if not pairs:
        return 1

    n_compare = len(pairs) if args.max_pairs <= 0 else min(args.max_pairs, len(pairs))

    n_pass = 0
    n_fail = 0

    for i, (fbin, fnc4) in enumerate(pairs[:n_compare], start=1):
        issues, n_bad, n_obs = compare_one(fbin, fnc4, rtol=args.rtol, atol=args.atol)
        rel = fbin.relative_to(root)
        if len(issues) == 0:
            print(f'[{i}/{n_compare}] PASS {rel} (N_obs={n_obs})')
            n_pass += 1
        else:
            print(f'[{i}/{n_compare}] FAIL {rel} (N_obs={n_obs}, n_bad={n_bad})')
            for msg in issues:
                print(f'  - {msg}')
            n_fail += 1

    print(f'summary: pass={n_pass}, fail={n_fail}, compared={n_compare}')
    return 0 if n_fail == 0 else 2


if __name__ == '__main__':
    raise SystemExit(main())
