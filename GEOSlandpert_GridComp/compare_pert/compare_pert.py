#!/usr/bin/env python3
"""
compare_pert.py

Read latlon_pert.bin and sphere_pert.bin (both raw Fortran binary streams
of real(4), shape 720 x 360, column-major = lon-fastest), and produce:

  1. Summary statistics table
  2. Zonal-mean profiles (mean and std vs latitude)
  3. Empirical covariance as a function of angular separation for
     a sample of reference points
  4. Side-by-side global maps (lat-lon pert and sphere pert)
  5. Difference map

Output:  compare_pert.png   (single figure with 5 panels)
         compare_pert.txt   (plain-text statistics)

Usage:
  python3 compare_pert.py

Requires: numpy, matplotlib (cartopy optional for nicer maps)
"""

import numpy as np
import os
import sys

# ── parameters matching the Fortran programs ─────────────────────────────────
N_LON  = 720
N_LAT  = 360
XCORR  = 10.0   # [degrees]

# ── load binary fields ────────────────────────────────────────────────────────
def load_field(fname):
    """Load a flat Fortran stream binary of N_LON*N_LAT real(4) values."""
    data = np.fromfile(fname, dtype=np.float32)
    if data.size != N_LON * N_LAT:
        raise ValueError(f"{fname}: expected {N_LON*N_LAT} values, got {data.size}")
    # Fortran column-major: first index (lon) varies fastest
    return data.reshape(N_LON, N_LAT, order='F')   # shape (lon, lat)

for fname in ('latlon_pert.bin', 'sphere_pert.bin'):
    if not os.path.exists(fname):
        sys.exit(f"ERROR: {fname} not found.  Run both MPI programs first.")

ll = load_field('latlon_pert.bin')   # (720, 360)
sp = load_field('sphere_pert.bin')   # (720, 360)

lons = np.linspace(-180 + 0.5, 180 - 0.5, N_LON)
lats = np.linspace(-90  + 0.5, 90  - 0.5, N_LAT)
LON2D, LAT2D = np.meshgrid(lons, lats, indexing='ij')   # (360,180)

# ── statistics ────────────────────────────────────────────────────────────────
def stats(f, name):
    return dict(name=name,
                mean=float(f.mean()), std=float(f.std()),
                min=float(f.min()),   max=float(f.max()))

st = {k: stats(f, k) for k, f in [('latlon', ll), ('sphere', sp)]}

lines = []
lines.append(f"{'':12s}  {'mean':>10s}  {'std':>10s}  {'min':>10s}  {'max':>10s}")
lines.append('-' * 60)
for k, s in st.items():
    lines.append(f"{s['name']:12s}  {s['mean']:10.5f}  {s['std']:10.5f}"
                 f"  {s['min']:10.5f}  {s['max']:10.5f}")
diff = ll - sp
d_st = stats(diff, 'diff (ll-sp)')
lines.append(f"\n{d_st['name']:12s}  {d_st['mean']:10.5f}  {d_st['std']:10.5f}"
             f"  {d_st['min']:10.5f}  {d_st['max']:10.5f}")

with open('compare_pert.txt', 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
print('\n'.join(lines))

# ── empirical zonal-mean and zonal-std vs latitude ───────────────────────────
ll_zmean = ll.mean(axis=0)   # (180,)
sp_zmean = sp.mean(axis=0)
ll_zstd  = ll.std(axis=0)
sp_zstd  = sp.std(axis=0)

# ── empirical covariance vs angular distance ──────────────────────────────────
# Sample a few reference points; compute covariance along a meridional transect.
def empirical_cov_meridional(field, ref_lat_idx, n_lons=36):
    """
    Average the covariance between a reference row (ref_lat_idx) and all
    other rows, using n_lons equally spaced longitudes.
    Returns (delta_lat_deg, cov) arrays.
    """
    lon_step = N_LON // n_lons
    cov = np.zeros(N_LAT)
    ref = field[::lon_step, ref_lat_idx]   # (n_lons,)
    for j in range(N_LAT):
        cov[j] = np.mean(ref * field[::lon_step, j])
    # normalise by cov[0]
    if cov[ref_lat_idx] > 0:
        cov /= cov[ref_lat_idx]
    dlat = lats - lats[ref_lat_idx]
    return dlat, cov

ref_idx = N_LAT // 2   # equator
dlat_ll, cov_ll = empirical_cov_meridional(ll, ref_idx)
dlat_sp, cov_sp = empirical_cov_meridional(sp, ref_idx)

# ── plotting ──────────────────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
except ImportError:
    sys.exit("matplotlib not available; statistics written to compare_pert.txt")

fig = plt.figure(figsize=(18, 14))
fig.suptitle(f'Perturbation comparison: flat lat-lon FFT vs sphere SMERFS\n'
             f'Grid {N_LON}x{N_LAT}, xcorr = ycorr = {XCORR}°',
             fontsize=13, fontweight='bold')

cmap   = 'RdBu_r'
vmax   = max(abs(ll).max(), abs(sp).max()) * 0.8
levels = np.linspace(-vmax, vmax, 41)

# -- Panel 1: lat-lon FFT map ------------------------------------------------
ax1 = fig.add_subplot(3, 3, 1)
c1 = ax1.contourf(lons, lats, ll.T, levels=levels, cmap=cmap, extend='both')
fig.colorbar(c1, ax=ax1, shrink=0.8)
ax1.set_title('Flat lat-lon FFT')
ax1.set_xlabel('Longitude [°]')
ax1.set_ylabel('Latitude [°]')

# -- Panel 2: sphere SMERFS map -----------------------------------------------
ax2 = fig.add_subplot(3, 3, 2)
c2 = ax2.contourf(lons, lats, sp.T, levels=levels, cmap=cmap, extend='both')
fig.colorbar(c2, ax=ax2, shrink=0.8)
ax2.set_title('Sphere SMERFS')
ax2.set_xlabel('Longitude [°]')

# -- Panel 3: difference -------------------------------------------------------
ax3 = fig.add_subplot(3, 3, 3)
dlev = np.linspace(-vmax*0.5, vmax*0.5, 41)
c3 = ax3.contourf(lons, lats, (ll - sp).T, levels=dlev, cmap='bwr', extend='both')
fig.colorbar(c3, ax=ax3, shrink=0.8)
ax3.set_title('Difference (FFT − Sphere)')
ax3.set_xlabel('Longitude [°]')

# -- Panel 4: zonal-mean profile ----------------------------------------------
ax4 = fig.add_subplot(3, 3, 4)
ax4.plot(ll_zmean, lats, 'b-',  label='FFT zonal mean')
ax4.plot(sp_zmean, lats, 'r--', label='Sphere zonal mean')
ax4.axvline(0, color='k', lw=0.5)
ax4.set_xlabel('Zonal mean')
ax4.set_ylabel('Latitude [°]')
ax4.set_title('Zonal mean vs latitude')
ax4.legend(fontsize=8)
ax4.grid(True, lw=0.3)

# -- Panel 5: zonal-std profile -----------------------------------------------
ax5 = fig.add_subplot(3, 3, 5)
ax5.plot(ll_zstd, lats, 'b-',  label='FFT zonal std')
ax5.plot(sp_zstd, lats, 'r--', label='Sphere zonal std')
ax5.set_xlabel('Zonal std')
ax5.set_ylabel('Latitude [°]')
ax5.set_title('Zonal std vs latitude')
ax5.legend(fontsize=8)
ax5.grid(True, lw=0.3)

# -- Panel 6: empirical covariance vs angular distance ------------------------
ax6 = fig.add_subplot(3, 3, 6)
efold_line = np.exp(-np.abs(dlat_ll) / XCORR)
ax6.plot(dlat_ll, cov_ll, 'b-',  label='FFT (meridional)')
ax6.plot(dlat_sp, cov_sp, 'r--', label='Sphere (meridional)')
ax6.plot(dlat_ll, efold_line, 'k:', lw=1, label=f'exp(-|Δlat|/{XCORR}°)')
ax6.axhline(np.exp(-1), color='gray', lw=0.7, ls='--', label='e⁻¹ level')
ax6.axhline(0, color='k', lw=0.5)
ax6.set_xlabel('Δ latitude [°]')
ax6.set_ylabel('Normalised covariance')
ax6.set_title(f'Meridional covariance (ref: equator)')
ax6.legend(fontsize=8)
ax6.set_xlim(-90, 90)
ax6.grid(True, lw=0.3)

# -- Panel 7: power spectrum (variance per latitude band) --------------------
ax7 = fig.add_subplot(3, 3, 7)
ax7.plot(ll_zstd**2, lats, 'b-',  label='FFT variance')
ax7.plot(sp_zstd**2, lats, 'r--', label='Sphere variance')
ax7.set_xlabel('Zonal variance')
ax7.set_ylabel('Latitude [°]')
ax7.set_title('Variance vs latitude (isotropy test)')
ax7.legend(fontsize=8)
ax7.grid(True, lw=0.3)

# -- Panel 8: histogram comparison -------------------------------------------
ax8 = fig.add_subplot(3, 3, 8)
bins = np.linspace(-4, 4, 60)
ax8.hist(ll.ravel(), bins=bins, alpha=0.6, color='blue',
         density=True, label='FFT')
ax8.hist(sp.ravel(), bins=bins, alpha=0.6, color='red',
         density=True, label='Sphere')
xg = np.linspace(-4, 4, 200)
from math import pi, exp
ax8.plot(xg, [1/((2*pi)**0.5)*exp(-x**2/2) for x in xg],
         'k-', lw=1.5, label='N(0,1)')
ax8.set_xlabel('Value')
ax8.set_ylabel('PDF')
ax8.set_title('Marginal distribution')
ax8.legend(fontsize=8)
ax8.grid(True, lw=0.3)

# -- Panel 9: scatter plot ---------------------------------------------------
ax9 = fig.add_subplot(3, 3, 9)
skip = 4   # downsample for readability
ax9.scatter(ll[::skip, ::skip].ravel(), sp[::skip, ::skip].ravel(),
            s=1, alpha=0.2, color='purple')
lim = vmax
ax9.plot([-lim, lim], [-lim, lim], 'k-', lw=0.8)
ax9.set_xlabel('FFT value')
ax9.set_ylabel('Sphere value')
ax9.set_title('FFT vs Sphere (scatter, decimated)')
ax9.set_xlim(-lim, lim)
ax9.set_ylim(-lim, lim)
ax9.set_aspect('equal')
ax9.grid(True, lw=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('compare_pert.png', dpi=130, bbox_inches='tight')
print('\nFigure saved to compare_pert.png')
