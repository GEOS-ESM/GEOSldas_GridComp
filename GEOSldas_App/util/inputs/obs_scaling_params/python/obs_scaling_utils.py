
# Helpers for the moving-window loop over the climatological year in the GEOSldas
# observation scaling scripts.  The GEOSldas output files are read with
# ../../shared/python/read_GEOSldas.py.
#
# amfox - Oct 2026

import os

import numpy as np

from datetime import datetime

nodata = -9999.

# ----------------------------------------------------------------------------
#
# helpers for the loop over the climatological year

# Days per month.  The climatological year is always the non-leap year 2014.

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

def get_ObsFcstAna_fname(inpath, exp_run, year, month, day, hour, minute):

    return (inpath + '/ana/ens_avg/Y' + f'{year:04d}' + '/M' + f'{month:02d}' + '/' + exp_run +
            '.ens_avg.ldas_ObsFcstAna.' + f'{year:04d}{month:02d}{day:02d}_{hour:02d}{minute:02d}z.nc4')

def any_ObsFcstAna_file(inpath, exp_run, run_months, first_year, last_year, cycles):

    # True if an ObsFcstAna file exists for run_months, years and cycles (seconds of day);
    # stops at the first file found

    for imonth, month in enumerate(run_months):
        for day in range(1, days_in_month[month-1] + 1):
            for seconds in cycles:
                for year in range(first_year[imonth], last_year[imonth] + 1):
                    fname = get_ObsFcstAna_fname(inpath, exp_run, year, month, day,
                                                 seconds // 3600, (seconds % 3600) // 60)
                    if os.path.exists(fname):
                        return True
    return False

def get_year_bounds(run_months, start_year, start_month, end_year, end_month):

    # first and last year to read for each entry of run_months, so that only
    # months between (start_year, start_month) and (end_year, end_month) are used

    first_year = []
    last_year  = []
    for month in run_months:
        first_year.append(start_year + 1 if month < start_month else start_year)
        last_year.append( end_year   - 1 if month > end_month   else end_year  )

    return first_year, last_year

def get_doy_pentad_range(run_months, first_year, last_year):

    # first/last DOY and pentad covered by run_months, for the output file names

    run_months = np.asarray(run_months)
    first_year = np.asarray(first_year)
    last_year  = np.asarray(last_year)

    first_month = int(run_months[first_year == first_year.min()].min())
    last_month  = int(run_months[last_year  == last_year.max() ].max())

    cumulative_days = np.concatenate(([0], np.cumsum(days_in_month)))

    doy_start    = int(cumulative_days[first_month-1] + 1)
    doy_end      = int(cumulative_days[last_month])
    pentad_start = int(np.ceil( doy_start / 5))
    pentad_end   = int(np.floor(doy_end   / 5))

    return doy_start, doy_end, pentad_start, pentad_end

def days_since_1950(times):

    return np.array([(t - datetime(1950, 1, 1)).total_seconds() / 86400. for t in times])

# ====================== EOF =========================================================
