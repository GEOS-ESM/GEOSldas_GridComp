#!/usr/bin/env python3
"""
pole_check.py -- row-free xcorr estimate for the sphere field in several regions.

Samples points uniformly by area inside a region, correlates all pairs within
2.5*xcorr, bins the products by great-circle distance, and fits
C(d) = exp(-d^2 / (2 sigma^2)). The field is normalised with the region's own
mean/std, so there are no latitude rows and no per-row means. Used to show
that the polar deviation of the row-based estimate (xcorr_stats.py) is an
estimator artifact, not a generator error.

Usage:  python3 pole_check.py [nlon nlat [xcorr ...]]     (needs scipy)
"""
import sys

import numpy as np
from scipy.spatial import cKDTree

NLON = int(sys.argv[1]) if len(sys.argv) > 2 else 2880
NLAT = int(sys.argv[2]) if len(sys.argv) > 2 else 1440
XCORRS = sys.argv[3:] if len(sys.argv) > 3 else ['0.2', '0.4', '0.8', '1.0']
DLON, DLAT = 360.0 / NLON, 180.0 / NLAT
LON = -180 + DLON * (np.arange(NLON) + 0.5)
LAT = -90 + DLAT * (np.arange(NLAT) + 0.5)
LOd, LAd = np.meshgrid(LON, LAT, indexing='ij')
LOr, LAr = np.deg2rad(LOd).ravel(order='F'), np.deg2rad(LAd).ravel(order='F')
rng = np.random.default_rng(0)

REGIONS = [
    ('eq(5S-5N,60lon)', (np.abs(LAd) < 5) & (np.abs(LOd) < 30)),
    ('45N(40-50N,80lon)', (np.abs(LAd - 45) < 5) & (np.abs(LOd) < 40)),
    ('cap>80N', LAd > 80),
    ('cap<80S', LAd < -80),
    ('cap>85N', LAd > 85),
]


def fit_sigma(f, mask, x, n=40000, nbins=25):
    idx = np.flatnonzero(mask.ravel(order='F'))
    w = np.cos(LAr[idx]); w /= w.sum()                       # area weights
    s = rng.choice(idx, min(n, len(idx) // 2), replace=False, p=w)
    lo, la = LOr[s], LAr[s]
    v = f.ravel(order='F')[s]
    v = (v - v.mean()) / v.std()
    P = np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
    rmax = 2.5 * x
    pr = cKDTree(P).query_pairs(2 * np.sin(np.deg2rad(rmax) / 2), output_type='ndarray')
    dist = np.rad2deg(2 * np.arcsin(np.linalg.norm(P[pr[:, 0]] - P[pr[:, 1]], axis=1) / 2))
    prod = v[pr[:, 0]] * v[pr[:, 1]]
    edges = np.linspace(0, rmax, nbins + 1)
    ib = np.digitize(dist, edges) - 1
    cnt = np.bincount(ib, minlength=nbins)[:nbins]
    c = np.bincount(ib, prod, minlength=nbins)[:nbins] / np.where(cnt > 0, cnt, np.nan)
    r = 0.5 * (edges[:-1] + edges[1:])
    sig = np.linspace(0.3 * x, 2 * x, 400)
    err = [np.nansum((c - np.exp(-r ** 2 / (2 * q * q))) ** 2) for q in sig]
    return sig[np.argmin(err)]


def main():
    print('sphere field: xcorr fitted from point pairs vs great-circle distance '
          '(no rows, no per-row mean)')
    print(f"{'xcorr':>5} " + ''.join(f'{name:>19s}' for name, _ in REGIONS))
    for xc in XCORRS:
        x = float(xc)
        f = np.fromfile(f'x{xc}/sphere_pert.bin', np.float32).reshape(NLON, NLAT, order='F')
        print(f'{xc:>5} ' + ''.join(f'{fit_sigma(f, m, x):19.3f}' for _, m in REGIONS))


if __name__ == '__main__':
    main()
