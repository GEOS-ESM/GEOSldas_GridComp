#!/usr/bin/env python3

# Jointly-matched OL legacy-vs-H121 ASCAT comparison, single experiment
# (OLv7_M36_MULTI_type_13_H121), no DA run involved.
#
# For each platform (Metop-A/B/C), a tile+cycle only contributes to the legacy
# species' sums if the H121 species also has a valid obs on that same tile in that
# same cycle, and vice versa (postproc_ObsFcstAna_selfxmask.py). This gives O_mean,
# F_mean, OmF_mean, OmF_stdv, N_data etc. for all 6 ASCAT species on a common,
# jointly-matched sample, so legacy vs H121 can be differenced directly for any
# metric -- unlike OL_vs_legacyobs / OL_vs_H121obs, which are each cross-masked
# independently against their own DA run and therefore draw on different
# day-subsets of the same OL trajectory (see raw-sample audit, July 2026: strict
# tile+cycle matched F_mean delta ~= 0 while the independently-sampled products
# show a spurious delta_F on par with delta_O).
#
# Usage (on Discover):
#   module load python/GEOSpyD
#   ./run_h121_legacy_xmask.py

import sys;  sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")
import os

from datetime import datetime
from netCDF4 import Dataset
import numpy as np

from run_h121_6yr import load_exp
from postproc_ObsFcstAna_selfxmask import postproc_ObsFcstAna_SelfXMask

DOMAIN = 'SMAP_EASEv2_M36_GLOBAL'
SUM_ROOT = '/gpfsm/dnb06/projects/p284/hsaf_cdr_test/omf_compare_sums/'
EXPTAG = 'OL_legacy_h121_xmask'

# obsparam 'species' values: legacy Metop-A/B/C (5,6,7) <-> H121 Metop-A/B/C (8,9,10)
SPECIES_PAIRS = {5: 8, 6: 9, 7: 10}

start_time = datetime(2015, 4, 1)
end_time = datetime(2021, 4, 1)   # exclusive -> matches OL_vs_H121obs / OL_vs_legacyobs span


def write_monthly_stats_nc4(file_path, stats):
    n_month, n_spec = stats['N_data'].shape
    with Dataset(file_path, 'w', format='NETCDF4') as nc:
        nc.createDimension('month', n_month)
        nc.createDimension('species', n_spec)

        yyyymm = nc.createVariable('yyyymm', 'i4', ('month',))
        yyyymm[:] = np.array([int(d) for d in stats['date_vec']])

        for key, value in stats.items():
            if key == 'date_vec':
                continue
            var = nc.createVariable(key, 'f4', ('month', 'species'),
                                     fill_value=-9999.0, zlib=True, complevel=4)
            var[:, :] = value


def main():
    sum_path = SUM_ROOT + EXPTAG + '/'
    os.makedirs(sum_path, exist_ok=True)

    ol_exp = load_exp('OLv7_M36_MULTI_type_13_H121', exptag=EXPTAG)

    postproc = postproc_ObsFcstAna_SelfXMask(ol_exp, SPECIES_PAIRS, start_time, end_time,
                                              sum_path=sum_path)
    postproc.save_monthly_sums()

    temporal_stats = postproc.calc_temporal_stats_from_sums(
        write_to_nc=True, fout_stats=sum_path + f'{EXPTAG}_temporal_stats.nc4')
    print(f'wrote {sum_path}{EXPTAG}_temporal_stats.nc4')

    monthly_stats = postproc.calc_spatial_stats_from_sums()
    write_monthly_stats_nc4(sum_path + f'{EXPTAG}_monthly_stats.nc4', monthly_stats)
    print(f'wrote {sum_path}{EXPTAG}_monthly_stats.nc4')


if __name__ == '__main__':
    main()

# ====================== EOF =========================================================
