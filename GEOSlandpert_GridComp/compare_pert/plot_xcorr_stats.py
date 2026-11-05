#!/usr/bin/env python3
"""
plot_xcorr_stats.py -- visualize xcorr recovered from the perturbations
(produces xcorr_stats.png).

Uses the row-based estimator from xcorr_stats.py; same inputs (x<xcorr>/ dirs).
  row 1 : east-west / north-south xcorr vs latitude, per target xcorr
  row 2 : band mean +- std of estimate / target
  row 3 : global summary vs xcorr, anisotropy ratio (E-W / N-S) for rect and
          sph, and a summary table

Usage:  python3 plot_xcorr_stats.py [nlon nlat [xcorr ...]]
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from xcorr_stats import LATS, XCORRS, NLON, NLAT, load, estimate

GENS = [('rect', 'latlon_pert.bin', 'tab:blue'), ('sph', 'sphere_pert.bin', 'tab:red')]
BANDS = [('|lat|<30', 0, 30), ('30-60', 30, 60), ('60-80', 60, 80), ('80-90', 80, 90)]


def smooth(v, n=17):
    k = np.ones(n)
    ok = ~np.isnan(v)
    return np.convolve(np.where(ok, v, 0), k, 'same') / np.maximum(np.convolve(ok, k, 'same'), 1)


def main():
    R = {(xc, g): estimate(load(f'x{xc}/{fn}'), float(xc)) for xc in XCORRS for g, fn, _ in GENS}
    nx = len(XCORRS)
    fig, axs = plt.subplots(3, max(nx, 4), figsize=(5 * max(nx, 4), 14))
    fig.suptitle(f'xcorr estimated from perturbations ({NLON}x{NLAT})  -- '
                 'rect = rectangle_random_fields, sph = sphere_random_fields',
                 fontsize=13, fontweight='bold')

    # row 1: estimates vs latitude
    for i, xc in enumerate(XCORRS):
        ax, x = axs[0, i], float(xc)
        for g, _, c in GENS:
            z, m = R[xc, g]
            ax.plot(LATS, smooth(z), c=c, label=f'{g} east-west')
            ax.plot(LATS, smooth(m), c=c, ls='--', label=f'{g} north-south')
        ax.plot(LATS, x * np.cos(np.deg2rad(LATS)), 'k:', lw=1, label='xcorr·cos(lat)')
        ax.axhline(x, c='k', lw=1.5, alpha=.4, label='target')
        ax.set_ylim(0, 1.6 * x); ax.set_xlim(-90, 90)
        ax.set_title(f'target xcorr = {xc}°'); ax.set_xlabel('lat')
        ax.set_ylabel('estimated xcorr [great-circle deg]'); ax.grid(lw=.3)
        if i == 0:
            ax.legend(fontsize=8)

    # row 2: band statistics
    w, xb = 0.2, np.arange(len(BANDS))
    for i, xc in enumerate(XCORRS):
        ax, x = axs[1, i], float(xc)
        for o, (g, _, c) in enumerate(GENS):
            for d_, (lab, which, hatch) in enumerate([('E-W', 0, None), ('N-S', 1, '//')]):
                v = R[xc, g][which]
                mu, sd = [], []
                for _, lo, hi in BANDS:
                    k = (np.abs(LATS) >= lo) & (np.abs(LATS) < hi)
                    q = v[k] / x
                    q = q[~np.isnan(q)]
                    mu.append(q.mean()); sd.append(q.std())
                ax.bar(xb + (2 * o + d_ - 1.5) * w, mu, w, yerr=sd, color=c,
                       alpha=.85 if which == 0 else .4, hatch=hatch, capsize=2, label=f'{g} {lab}')
        ax.axhline(1, c='k', lw=1)
        ax.set_xticks(xb); ax.set_xticklabels([b[0] for b in BANDS]); ax.set_ylim(0, 1.7)
        ax.set_ylabel('estimate / target (mean ± std over rows)')
        ax.set_title(f'band statistics, xcorr = {xc}°'); ax.grid(lw=.3, axis='y')
        if i == 0:
            ax.legend(fontsize=8, ncol=2)

    # row 3a: global summary
    ax = axs[2, 0]
    xs = [float(x) for x in XCORRS]
    for g, _, c in GENS:
        for which, ls in [(0, '-'), (1, '--')]:
            mu, p5, p95 = [], [], []
            for xc in XCORRS:
                v = R[xc, g][which][np.abs(LATS) < 80] / float(xc)
                v = v[~np.isnan(v)]
                mu.append(v.mean()); p5.append(np.percentile(v, 5)); p95.append(np.percentile(v, 95))
            ax.plot(xs, mu, 'o' + ls, c=c, label=f"{g} {'E-W' if which == 0 else 'N-S'}")
            ax.fill_between(xs, p5, p95, color=c, alpha=.1)
    ax.axhline(1, c='k', lw=1); ax.set_xlabel('target xcorr [deg]'); ax.set_ylabel('estimate / target')
    ax.set_title('global |lat|<80: mean, 5–95% shaded'); ax.legend(fontsize=8); ax.grid(lw=.3)

    # row 3b,c: anisotropy ratio
    for j, (g, _, _) in enumerate(GENS):
        ax = axs[2, 1 + j]
        for xc, cc in zip(XCORRS, plt.cm.viridis(np.linspace(0, .9, nx))):
            z, m = R[xc, g]
            ax.plot(LATS, smooth(z / m, 33), c=cc, label=f'xcorr={xc}')
        ax.plot(LATS, np.cos(np.deg2rad(LATS)), 'k:', label='cos(lat)')
        ax.axhline(1, c='k', lw=1.5, alpha=.4)
        ax.set_ylim(0, 1.5); ax.set_xlim(-90, 90); ax.set_xlabel('lat')
        ax.set_ylabel('E-W / N-S xcorr (anisotropy)')
        ax.set_title(f'{g}: anisotropy ratio (target 1)'); ax.legend(fontsize=8); ax.grid(lw=.3)

    # row 3d: table
    ax = axs[2, 3]; ax.axis('off')
    rows = []
    for xc in XCORRS:
        for g, _, _ in GENS:
            z, m = R[xc, g]
            k = np.abs(LATS) < 80
            rows.append([xc, g, f'{np.nanmean(z[k]):.3f}±{np.nanstd(z[k]):.3f}',
                         f'{np.nanmean(m[k]):.3f}±{np.nanstd(m[k]):.3f}', f'{np.nanmedian(z[k] / m[k]):.2f}'])
    t = ax.table(cellText=rows, colLabels=['target', 'gen', 'E-W xcorr', 'N-S xcorr', 'E-W/N-S'],
                 loc='center', cellLoc='center')
    t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.5)
    for r in range(1, len(rows) + 1):
        t[r, 1].set_facecolor('#cfe2f3' if rows[r - 1][1] == 'rect' else '#f4cccc')
    ax.set_title('global |lat|<80 (deg)')
    for ax in axs[2, 4:]:
        ax.axis('off')

    plt.tight_layout(rect=[0, 0, 1, .96])
    plt.savefig('xcorr_stats.png', dpi=110)
    print('Figure saved to xcorr_stats.png')


if __name__ == '__main__':
    main()
