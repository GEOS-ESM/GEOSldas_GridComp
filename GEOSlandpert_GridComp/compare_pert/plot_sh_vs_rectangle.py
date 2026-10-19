#!/usr/bin/env python3
"""
plot_sh_vs_rectangle.py

Visualise the perturbation fields written by sh_vs_rectangle.x:

  sh_pert.bin            sphere_harmonics_random_fields
  rect_pert.bin          rectangle_random_fields (coarsen = .false.)
  rect_coarsen_pert.bin  rectangle_random_fields (coarsen = .true.), optional

File format: int32 nx, int32 ny, then real(4) f(nx, ny) (Fortran order,
lon fastest, lat from the south pole), cell-centred global grid.

Figure rows:
  1. global maps
  2. 20x20 deg zooms at the equator and at 60N, drawn with true km aspect
     (x scaled by cos(lat)) so the flat-FFT anisotropy is visible
  3. correlation vs great-circle distance (km): E-W at equator, E-W at
     60N, N-S at equator; Gaussian exp(-d^2/2s^2), s = xcorr, for reference

Usage:  python3 plot_sh_vs_rectangle.py [xcorr_deg] [output.png]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

XCORR = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
OUT = sys.argv[2] if len(sys.argv) > 2 else 'sh_vs_rectangle.png'
R_KM = 6371.0
KM_PER_DEG = np.pi * R_KM / 180.0


def load(fname):
    with open(fname, 'rb') as fh:
        nx, ny = np.fromfile(fh, dtype=np.int32, count=2)
        f = np.fromfile(fh, dtype=np.float32, count=nx * ny)
    return f.reshape(nx, ny, order='F').T          # (lat, lon), south first


cases = [('SH (sphere harmonics)', 'sh_pert.bin'),
         ('RECT full (flat FFT)', 'rect_pert.bin'),
         ('RECT coarsen', 'rect_coarsen_pert.bin')]
cases = [(t, load(f)) for t, f in cases if os.path.exists(f)]
if len(cases) < 2:
    sys.exit('Need sh_pert.bin and rect_pert.bin; run sh_vs_rectangle.x first.')


def axes_of(f):
    ny, nx = f.shape
    lon = -180 + (np.arange(nx) + 0.5) * 360.0 / nx
    lat = -90 + (np.arange(ny) + 0.5) * 180.0 / ny
    return lon, lat


def corr_vs_dist(f, lat0, along_lat, maxdeg):
    """Lagged correlation in a +-5 deg band around lat0."""
    ny, nx = f.shape
    dy, dx = 180.0 / ny, 360.0 / nx
    j0 = int(round((lat0 + 90) / dy - 0.5))
    jb = max(1, int(round(5 / dy)))
    band = f[max(0, j0 - jb):j0 + jb + 1, :].astype(np.float64)
    nlag = int(maxdeg / (dy if along_lat else dx)) + 1
    d, c = [], []
    for lag in range(nlag):
        if along_lat:
            a = f[max(0, j0 - jb):j0 + jb + 1, :].astype(np.float64)
            b = f[max(0, j0 - jb) + lag:j0 + jb + 1 + lag, :].astype(np.float64)
            n = min(a.shape[0], b.shape[0]); a, b = a[:n], b[:n]
            d.append(lag * dy * KM_PER_DEG)
        else:
            a, b = band, np.roll(band, -lag, axis=1)
            d.append(lag * dx * KM_PER_DEG * np.cos(np.radians(lat0)))
        c.append((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))
    return np.array(d), np.array(c)


ncol = len(cases)
fig = plt.figure(figsize=(6 * max(ncol, 3), 16))
fig.suptitle(f'Perturbation fields, xcorr = {XCORR} deg  '
             f'(first field of first draw)', fontsize=14, fontweight='bold')
vmax = 3.0

# --- row 1: global maps -----------------------------------------------------
for i, (title, f) in enumerate(cases):
    lon, lat = axes_of(f)
    ax = fig.add_subplot(4, ncol, i + 1)
    im = ax.imshow(f, origin='lower', extent=[-180, 180, -90, 90],
                   cmap='RdBu_r', vmin=-vmax, vmax=vmax, interpolation='nearest')
    ax.set_title(f'{title}  {f.shape[1]}x{f.shape[0]}')
    ax.set_xlabel('lon'); ax.set_ylabel('lat')
    fig.colorbar(im, ax=ax, shrink=0.7)

# --- rows 2-3: zooms with true km aspect -------------------------------------
box = max(20.0, 15 * XCORR)
for r, lat0 in enumerate([0.0, 60.0]):
    for i, (title, f) in enumerate(cases):
        lon, lat = axes_of(f)
        jj = (lat >= lat0 - box / 2) & (lat <= lat0 + box / 2)
        ii = (lon >= -box / 2) & (lon <= box / 2)
        ax = fig.add_subplot(4, ncol, (r + 1) * ncol + i + 1)
        ax.imshow(f[np.ix_(jj, ii)], origin='lower',
                  extent=[-box / 2, box / 2, lat0 - box / 2, lat0 + box / 2],
                  cmap='RdBu_r', vmin=-vmax, vmax=vmax, interpolation='nearest',
                  aspect=1.0 / np.cos(np.radians(lat0)))
        ax.set_title(f'{title.split(" (")[0]} zoom @ {lat0:.0f}N (km-true aspect)')
        ax.set_xlabel('lon'); ax.set_ylabel('lat')

# --- row 4: correlation vs distance ------------------------------------------
maxdeg = 4 * XCORR
sig = XCORR * KM_PER_DEG
panels = [('E-W @ equator', 0.0, False), ('E-W @ 60N', 60.0, False),
          ('N-S @ equator', 0.0, True)]
styles = ['r-', 'b--', 'g:']
for p, (ptitle, lat0, along_lat) in enumerate(panels):
    ax = fig.add_subplot(4, 3, 10 + p)
    for (title, f), st in zip(cases, styles):
        d, c = corr_vs_dist(f, lat0, along_lat,
                            maxdeg / (np.cos(np.radians(lat0)) if not along_lat else 1))
        ax.plot(d, c, st, marker='.', ms=3, label=title.split(' (')[0])
    dd = np.linspace(0, 4 * sig, 200)
    ax.plot(dd, np.exp(-dd ** 2 / (2 * sig ** 2)), 'k-', lw=0.8,
            label=f'Gaussian, sigma={sig:.0f} km')
    ax.axhline(np.exp(-0.5), color='gray', lw=0.6, ls='--')
    ax.axvline(sig, color='gray', lw=0.6, ls='--')
    ax.set_xlim(0, 4 * sig); ax.set_ylim(-0.2, 1.02)
    ax.set_xlabel('great-circle distance [km]'); ax.set_ylabel('correlation')
    ax.set_title(ptitle); ax.grid(True, lw=0.3); ax.legend(fontsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(OUT, dpi=110)
print(f'Figure saved to {OUT}')
