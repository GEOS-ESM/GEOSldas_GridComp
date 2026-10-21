#!/usr/bin/env python3
"""
compare_pert.py -- visualize perturbations from compare_pert.x

Reads latlon_pert.bin (rectangle_random_fields, flat lat-lon FFT) and
sphere_pert.bin (sphere_random_fields, spherical harmonics), plus
compare_pert_grid.txt (nlon nlat xcorr) written by compare_pert.x.

Figure panels:
  1-2  global maps (equirectangular)
  3-4  polar stereographic-style views of the north pole (lon/lat -> x,y)
  5    zonal std vs latitude (isotropy / pole behaviour)
  6    empirical correlation vs great-circle distance at several latitudes
  7    marginal histogram vs N(0,1)
  8    zonal (east-west) correlation in km at 0, 45, 70 deg -> shows that the
       flat FFT keeps a fixed angular length, i.e. shrinking physical length
       toward the poles, while the sphere field is isotropic in km.

Both generators use xcorr as a Gaussian sigma: C(d) = exp(-d^2 / (2 xcorr^2)).

Usage:  python3 compare_pert.py [output-prefix]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PREFIX = sys.argv[1] if len(sys.argv) > 1 else 'compare_pert'
R_EARTH = 6371.0

if os.path.exists('compare_pert_grid.txt'):
    v = open('compare_pert_grid.txt').read().split()
    NLON, NLAT, XCORR = int(v[0]), int(v[1]), float(v[2])
else:
    NLON, NLAT, XCORR = 720, 360, 2.0


def load(fname):
    d = np.fromfile(fname, dtype=np.float32)
    if d.size != NLON * NLAT:
        sys.exit(f'{fname}: expected {NLON*NLAT} values, got {d.size}')
    return d.reshape(NLON, NLAT, order='F')  # (lon, lat)


fields = {'Rectangle (flat FFT)': load('latlon_pert.bin'),
          'Sphere (spherical harmonics)': load('sphere_pert.bin')}
colors = {'Rectangle (flat FFT)': 'tab:blue', 'Sphere (spherical harmonics)': 'tab:red'}

dlon, dlat = 360.0 / NLON, 180.0 / NLAT
lons = -180 + dlon * (np.arange(NLON) + 0.5)
lats = -90 + dlat * (np.arange(NLAT) + 0.5)

# ---------------------------------------------------------------- statistics
lines = [f'grid {NLON}x{NLAT}, xcorr={XCORR} deg (Gaussian sigma)',
         f"{'':30s}{'mean':>10s}{'std':>10s}{'min':>10s}{'max':>10s}"]
for k, f in fields.items():
    lines.append(f'{k:30s}{f.mean():10.4f}{f.std():10.4f}{f.min():10.4f}{f.max():10.4f}')
open(f'{PREFIX}.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))


# ------------------------------------------------- correlation diagnostics
def zonal_corr(f, j, maxlag):
    """East-west correlation along latitude row j (periodic), lags 0..maxlag."""
    row = f[:, j] - f[:, j].mean()
    var = (row * row).mean()
    return np.array([(row * np.roll(row, -k)).mean() / var for k in range(maxlag + 1)])


def meridional_corr(f, j0, maxlag):
    """North-south correlation from row j0 to rows j0..j0+maxlag (avg over lon)."""
    a = f[:, j0] - f[:, j0].mean()
    out = []
    for k in range(maxlag + 1):
        b = f[:, j0 + k] - f[:, j0 + k].mean()
        out.append((a * b).mean() / np.sqrt((a * a).mean() * (b * b).mean()))
    return np.array(out)


def rows_near(lat0, half_width=5.0):
    return np.where(np.abs(lats - lat0) <= half_width)[0]


maxdeg = 4 * XCORR
fig = plt.figure(figsize=(18, 15))
fig.suptitle(f'rectangle_random_fields vs sphere_random_fields   '
             f'grid {NLON}x{NLAT}, xcorr = {XCORR}°', fontsize=14, fontweight='bold')

vmax = 3.0
# maps
for i, (k, f) in enumerate(fields.items()):
    ax = fig.add_subplot(4, 2, i + 1)
    im = ax.pcolormesh(lons, lats, f.T, cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
    ax.set_title(k)
    ax.set_xlabel('lon'); ax.set_ylabel('lat')
    fig.colorbar(im, ax=ax, shrink=0.85)

# polar views (north of 60N), azimuthal equidistant projection
jn = lats >= 60
LON, LAT = np.meshgrid(np.deg2rad(lons), lats[jn], indexing='ij')
rr = 90 - LAT
X, Y = rr * np.sin(LON), -rr * np.cos(LON)
for i, (k, f) in enumerate(fields.items()):
    ax = fig.add_subplot(4, 2, 3 + i)
    im = ax.pcolormesh(X, Y, f[:, jn], cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='gouraud')
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f'{k}: north of 60°N (azimuthal equidistant)')
    fig.colorbar(im, ax=ax, shrink=0.85)

# zonal std vs lat
ax = fig.add_subplot(4, 2, 5)
for k, f in fields.items():
    ax.plot(f.std(axis=0), lats, color=colors[k], label=k)
ax.axvline(1, color='k', lw=0.6, ls=':')
ax.set_xlabel('zonal std'); ax.set_ylabel('lat'); ax.set_title('Std vs latitude (should be 1)')
ax.legend(fontsize=8); ax.grid(lw=0.3)

# histogram
ax = fig.add_subplot(4, 2, 6)
bins = np.linspace(-4.5, 4.5, 80)
for k, f in fields.items():
    ax.hist(f.ravel(), bins=bins, density=True, histtype='step', lw=1.5, color=colors[k], label=k)
xg = np.linspace(-4.5, 4.5, 200)
ax.plot(xg, np.exp(-xg**2 / 2) / np.sqrt(2 * np.pi), 'k--', lw=1, label='N(0,1)')
ax.set_title('Marginal distribution'); ax.legend(fontsize=8); ax.grid(lw=0.3)

# zonal correlation vs physical distance at several latitudes
ax = fig.add_subplot(4, 2, 7)
ls_lat = {0: '-', 45: '--', 70: ':'}
for k, f in fields.items():
    for lat0, ls in ls_lat.items():
        rows = rows_near(lat0)
        coslat = np.cos(np.deg2rad(lat0))
        maxlag = int(min(NLON // 2, np.ceil(maxdeg / (dlon * max(coslat, 0.05)))))
        c = np.mean([zonal_corr(f, j, maxlag) for j in rows], axis=0)
        dist_km = np.arange(maxlag + 1) * dlon * coslat * np.pi / 180 * R_EARTH
        ax.plot(dist_km, c, color=colors[k], ls=ls, label=f'{k.split()[0]} {lat0}°')
dkm = np.linspace(0, maxdeg * np.pi / 180 * R_EARTH, 200)
sig_km = XCORR * np.pi / 180 * R_EARTH
ax.plot(dkm, np.exp(-dkm**2 / (2 * sig_km**2)), 'k', lw=2, alpha=0.4, label='target Gaussian')
ax.set_xlim(0, dkm[-1]); ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('east-west distance [km]'); ax.set_ylabel('correlation')
ax.set_title('Zonal correlation vs physical distance (isotropy in km)')
ax.legend(fontsize=7, ncol=2); ax.grid(lw=0.3)

# meridional correlation (same in both if correct)
ax = fig.add_subplot(4, 2, 8)
maxlag = int(np.ceil(maxdeg / dlat))
for k, f in fields.items():
    for lat0, ls in ls_lat.items():
        rows = rows_near(lat0)
        rows = rows[rows + maxlag < NLAT]
        c = np.mean([meridional_corr(f, j, maxlag) for j in rows], axis=0)
        ax.plot(np.arange(maxlag + 1) * dlat, c, color=colors[k], ls=ls,
                label=f'{k.split()[0]} {lat0}°')
dd = np.linspace(0, maxdeg, 200)
ax.plot(dd, np.exp(-dd**2 / (2 * XCORR**2)), 'k', lw=2, alpha=0.4, label='target Gaussian')
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('north-south separation [deg]'); ax.set_ylabel('correlation')
ax.set_title('Meridional correlation'); ax.legend(fontsize=7, ncol=2); ax.grid(lw=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(f'{PREFIX}.png', dpi=120)
print(f'Figure saved to {PREFIX}.png')
