#!/usr/bin/env python3
"""
xcorr_stats.py -- recover xcorr from perturbations written by compare_pert.x

Expects one sub-directory per xcorr, named x<xcorr>/, each holding
latlon_pert.bin (rectangle_random_fields) and sphere_pert.bin
(sphere_random_fields), as produced by run_xcorr_stats.sh.

For every latitude row, xcorr is the distance at which the correlation drops
to exp(-1/2) (= sigma for C(d) = exp(-d^2 / (2 sigma^2))):
  east-west   : correlation along the row, converted to great-circle deg
  north-south : correlation between the row and the rows north of it,
                at the same longitude
Prints mean +- std [5th, 95th percentile] over rows in latitude bands.

NOTE: this row-based estimator is biased within a few degrees of the poles
(see XCORR_STATS.md and pole_check.py).

Usage:  python3 xcorr_stats.py [nlon nlat [xcorr ...]]
"""
import sys

import numpy as np

NLON = int(sys.argv[1]) if len(sys.argv) > 2 else 2880
NLAT = int(sys.argv[2]) if len(sys.argv) > 2 else 1440
XCORRS = sys.argv[3:] if len(sys.argv) > 3 else ['0.2', '0.4', '0.8', '1.0']
DLAT = 180.0 / NLAT
DLON = 360.0 / NLON
LATS = -90 + DLAT * (np.arange(NLAT) + 0.5)
GENS = [('rect', 'latlon_pert.bin'), ('sph', 'sphere_pert.bin')]
BANDS = [('|lat|<30', 0, 30), ('30-60', 30, 60), ('60-80', 60, 80),
         ('80-90', 80, 90), ('global<80', 0, 80)]
T = np.exp(-0.5)


def load(path):
    """Read a (nlon, nlat) float32 Fortran-ordered field."""
    return np.fromfile(path, np.float32).reshape(NLON, NLAT, order='F')


def crossing(c, step):
    """Distance where correlation c (sampled every `step`) first drops below exp(-1/2)."""
    i = np.argmax(c < T)
    if i == 0:
        return np.nan
    return step * (i - 1 + (c[i - 1] - T) / (c[i - 1] - c[i]))


def estimate(f, xcorr):
    """Per-row east-west and north-south xcorr estimates [great-circle deg]."""
    a = f - f.mean(axis=0)
    var = (a * a).mean(axis=0)
    F = np.fft.rfft(a, axis=0)
    ac = np.fft.irfft(F * np.conj(F), n=NLON, axis=0) / NLON / var  # zonal autocorr
    zon = np.full(NLAT, np.nan)
    mer = np.full(NLAT, np.nan)
    maxlag = int(5 * xcorr / DLAT) + 2
    for j in range(NLAT):
        zon[j] = crossing(ac[:NLON // 2, j], DLON * np.cos(np.deg2rad(LATS[j])))
        if j + maxlag < NLAT:
            c = (a[:, j:j + maxlag + 1] * a[:, [j]]).mean(0) / np.sqrt(var[j] * var[j:j + maxlag + 1])
            mer[j] = crossing(c, DLAT)
    return zon, mer


def summary(v):
    v = v[~np.isnan(v)]
    return f'{v.mean():5.3f}±{v.std():5.3f} [{np.percentile(v, 5):5.3f},{np.percentile(v, 95):5.3f}]'


def main():
    for xc in XCORRS:
        print(f'\n=== target xcorr = {xc} deg   (mean±std over rows [5th,95th pct], great-circle deg)')
        print(f"{'':5s}{'band':10s}{'zonal xcorr':>30s}{'meridional xcorr':>30s}{'zonal/merid':>14s}")
        for g, fn in GENS:
            f = load(f'x{xc}/{fn}')
            zon, mer = estimate(f, float(xc))
            kurt = ((f - f.mean()) ** 4).mean() / f.var() ** 2
            print(f'{g:5s}field std={f.std():.3f} kurtosis={kurt:.2f}')
            for b, lo, hi in BANDS:
                k = (np.abs(LATS) >= lo) & (np.abs(LATS) < hi)
                r = zon[k] / mer[k]
                print(f'{g:5s}{b:10s}{summary(zon[k]):>30s}{summary(mer[k]):>30s}{np.nanmedian(r):14.2f}')


if __name__ == '__main__':
    main()
