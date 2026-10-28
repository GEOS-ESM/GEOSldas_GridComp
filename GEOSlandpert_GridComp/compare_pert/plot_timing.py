#!/usr/bin/env python3
"""
plot_timing.py -- plot compare_pert_timing.csv produced by compare_pert.x.

Figure: init and mean draw time vs MPI ranks (one line per xcorr) and
vs xcorr (at the largest rank count), for both generators.

Usage:  python3 plot_timing.py [csv] [output-prefix]
"""
import csv
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else 'compare_pert_timing.csv'
PREFIX = sys.argv[2] if len(sys.argv) > 2 else 'compare_pert_timing'

# keep the last entry for each (generator, npes, xcorr, grid)
rows = {}
for r in csv.DictReader(open(CSV)):
    key = (r['generator'], int(r['npes']), float(r['xcorr']), r['nlon'], r['nlat'])
    rows[key] = (float(r['init_s']), float(r['draw_mean_s']))

gens = sorted({k[0] for k in rows})
xcorrs = sorted({k[2] for k in rows})
npes_all = sorted({k[1] for k in rows})
grid = sorted({(k[3], k[4]) for k in rows})
style = {'rectangle': ('tab:blue', 'o'), 'sphere': ('tab:red', 's')}

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
fig.suptitle(f'Performance: rectangle_random_fields vs sphere_random_fields, '
             f'grid {" / ".join("x".join(g) for g in grid)}', fontweight='bold')

for ax, idx, name in [(axes[0], 0, 'init (constructor)'), (axes[1], 1, 'draw (2 fields / call)')]:
    for g in gens:
        c, m = style.get(g, ('k', 'x'))
        for i, xc in enumerate(xcorrs):
            pts = sorted((k[1], v[idx]) for k, v in rows.items() if k[0] == g and k[2] == xc)
            if pts:
                ax.plot(*zip(*pts), color=c, marker=m, alpha=0.35 + 0.65 * (i + 1) / len(xcorrs),
                        label=f'{g} xcorr={xc:g}')
    ax.set_xscale('log', base=2); ax.set_yscale('log')
    ax.set_xlabel('MPI ranks'); ax.set_ylabel('seconds'); ax.set_title(name)
    ax.grid(which='both', lw=0.3); ax.legend(fontsize=7)

ax = axes[2]
npmax = npes_all[-1]
for g in gens:
    c, m = style.get(g, ('k', 'x'))
    for idx, ls, lab in [(0, '--', 'init'), (1, '-', 'draw')]:
        pts = sorted((k[2], v[idx]) for k, v in rows.items() if k[0] == g and k[1] == npmax)
        if pts:
            ax.plot(*zip(*pts), color=c, marker=m, ls=ls, label=f'{g} {lab}')
ax.set_yscale('log'); ax.set_xlabel('xcorr [deg]'); ax.set_ylabel('seconds')
ax.set_title(f'Time vs xcorr at {npmax} ranks'); ax.grid(which='both', lw=0.3); ax.legend(fontsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig(f'{PREFIX}.png', dpi=120)
print(f'Figure saved to {PREFIX}.png')

# text table
print(f"\n{'generator':10s}{'npes':>6s}{'xcorr':>8s}{'init[s]':>12s}{'draw[s]':>12s}")
for k in sorted(rows):
    print(f'{k[0]:10s}{k[1]:6d}{k[2]:8.2f}{rows[k][0]:12.4f}{rows[k][1]:12.4f}')
