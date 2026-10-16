
# Compute climatological mean and stdv of observations and model forecasts from
# GEOSldas ObsFcstAna files on an Earth-fixed global lat/lon grid, for use in the
# z-score scaling of observations before assimilation (e.g., ASCAT soil moisture).
#
# For each day of the climatological year, the statistics are computed from all
# obs/fcst pairs within a moving window of w_days days (across all years), and
# stored at the center of the window.  Output is written for each pentad and/or
# DOY, and/or as one file with all pentads.
#
# Python version of get_model_and_obs_clim_stats_latlon_grid.m.  The grid cell of
# each obs is computed with ceil() in single precision, as in the lookup of the
# scaling parameters in GEOSldas (scale_obs_sfmc_zscore).
#
# amfox - Oct 2026

import os
import sys;        sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../shared/python'))
import warnings

import numpy as np

from datetime          import datetime, timedelta
from netCDF4           import Dataset

from read_GEOSldas      import read_ObsFcstAna_nc4
from obs_scaling_utils import get_ObsFcstAna_fname, get_doy_pentad_range
from obs_scaling_utils import days_in_month, days_since_1950, count_ObsFcstAna_files

n_fields  = 7      # o_mean, o_std, m_mean, m_std, n_data, m_min, m_max
n_pentads = 73

def get_obs_scaling_latlon(exp_path, exp_run, domain, species, species_names, combine_species,
                           run_months, first_year, last_year, dt_assim, t0_assim,
                           resol, w_days, ndata_min, prefix, out_dir,
                           print_each_doy, print_each_pentad, print_all_pentads,
                           exclude_tiles=None):

    # exclude_tiles: optional array of tile numbers (1-based) whose obs are not used, e.g.
    #                tiles flagged by a QC that postdates the ObsFcstAna files.

    if combine_species:
        species_groups = [list(species)]
        species_labels = ['ALL']
    else:
        species_groups = [[s] for s in species]
        species_labels = list(species_names)

    n_groups = len(species_groups)

    inpath  = exp_path + '/' + exp_run + '/output/' + domain
    outpath = inpath + '/stats/' + out_dir
    os.makedirs(outpath, exist_ok=True)

    t0_assim = t0_assim % dt_assim
    cycles   = range(t0_assim, 86400, dt_assim)

    n_files = count_ObsFcstAna_files(inpath, exp_run, run_months, first_year, last_year, cycles)
    if n_files == 0:
        raise ValueError(f'no ObsFcstAna files for the requested period in {inpath}/ana/ens_avg')

    # global lat/lon grid; grid cell index = j*n_lon + i  (i: lon index, j: lat index)
    # (keep ll_lon, ll_lat, d_lon, d_lat as Python floats, so that the cell index of
    #  the float32 obs lat/lon is computed in single precision)

    ll_lon  = -180.
    ll_lat  =  -90.
    d_lon   = float(resol)
    d_lat   = float(resol)
    n_lon   = int(round(360. / d_lon))
    n_lat   = int(round(180. / d_lat))
    n_cells = n_lon * n_lat

    # sums for each grid cell and each day of the moving window

    shape  = (n_groups, n_cells, w_days)
    o_sum  = np.zeros(shape)
    o_sum2 = np.zeros(shape)
    m_sum  = np.zeros(shape)
    m_sum2 = np.zeros(shape)
    n_data = np.zeros(shape)
    m_min  = np.full( shape,  np.inf)
    m_max  = np.full( shape, -np.inf)

    if print_all_pentads:
        data_out = np.full((n_groups, n_fields, n_cells, n_pentads), np.nan)
    start_time_p = [datetime(2014, 1, 1)] * n_pentads
    end_time_p   = [datetime(2014, 1, 1)] * n_pentads

    # output file names

    d1, d2, p1, p2 = get_doy_pentad_range(run_months, first_year, last_year)

    fname_base_d = (outpath + '/' + prefix + f'{min(first_year)}_doy{d1}_{max(last_year)}_doy{d2}' +
                    f'_W_{w_days}d_Nmin_{ndata_min}')
    fname_base_p = (outpath + '/' + prefix + f'{min(first_year)}_p{p1}_{max(last_year)}_p{p2}' +
                    f'_W_{int(round(w_days / 5))}p_Nmin_{ndata_min}')

    grid = (n_lon, n_lat, ll_lon, ll_lat, d_lon, d_lat)

    n_read = 0                          # number of ObsFcstAna files read
    n_obs  = np.zeros(n_groups, int)    # number of obs used for each species group

    count = 0

    for imonth, month in enumerate(run_months):
        for day in range(1, days_in_month[month-1] + 1):

            print(f"Processing {month:02d}-{day:02d}")

            if count < w_days:
                count += 1
            idx = count - 1             # index of this day in the moving window

            for seconds in cycles:
                hour   = seconds // 3600
                minute = (seconds % 3600) // 60

                for year in range(first_year[imonth], last_year[imonth] + 1):

                    fname = get_ObsFcstAna_fname(inpath, exp_run, year, month, day, hour, minute)
                    if not os.path.exists(fname):
                        continue
                    ofa = read_ObsFcstAna_nc4(fname)
                    n_read += 1

                    # use only obs/fcst pairs (obs and fcst both valid)
                    valid = ~np.isnan(ofa['obs_fcst']) & ~np.isnan(ofa['obs_obs'])

                    for g, species_group in enumerate(species_groups):

                        mask = valid & np.isin(ofa['obs_species'], species_group)
                        if exclude_tiles is not None:
                            mask &= ~np.isin(ofa['obs_tilenum'], exclude_tiles)
                        if not np.any(mask):
                            continue

                        # grid cell of each obs (drop obs outside the grid)

                        i_idx  = np.ceil((ofa['obs_lon'][mask] - ll_lon) / d_lon).astype(int) - 1
                        j_idx  = np.ceil((ofa['obs_lat'][mask] - ll_lat) / d_lat).astype(int) - 1
                        inside = (i_idx >= 0) & (i_idx < n_lon) & (j_idx >= 0) & (j_idx < n_lat)

                        cells  = j_idx[inside] * n_lon + i_idx[inside]
                        obs    = ofa['obs_obs' ][mask][inside].astype(float)
                        fcst   = ofa['obs_fcst'][mask][inside].astype(float)

                        # add all obs/fcst to their grid cell (np.add.at adds repeated cells
                        # one by one)

                        n_obs[g] += len(cells)
                        np.add.at(    o_sum[ g, :, idx], cells, obs        )
                        np.add.at(    o_sum2[g, :, idx], cells, obs * obs  )
                        np.add.at(    n_data[g, :, idx], cells, 1.         )

                        np.add.at(    m_sum[ g, :, idx], cells, fcst       )
                        np.add.at(    m_sum2[g, :, idx], cells, fcst * fcst)
                        np.minimum.at(m_min[ g, :, idx], cells, fcst       )
                        np.maximum.at(m_max[ g, :, idx], cells, fcst       )

            if count < w_days:
                continue

            # The window is full: compute the statistics, store them at the center of the
            # window, and shift the window by one day.

            end_time   = datetime(2014, month, day, hour, minute)
            start_time = end_time - timedelta(seconds=w_days*86400)
            doy        = (end_time - timedelta(seconds=int(w_days*86400/2.))).timetuple().tm_yday
            pentad     = (doy + 2) // 5

            is_pentad_center = (doy + 2) % 5 == 0

            for g in range(n_groups):

                if not (print_each_doy or (is_pentad_center and (print_each_pentad or print_all_pentads))):
                    continue

                data2d = get_window_stats(o_sum[g], o_sum2[g], m_sum[g], m_sum2[g], n_data[g],
                                          m_min[g], m_max[g], ndata_min)

                if print_each_doy:
                    fname_out = fname_base_d + '_sp_' + species_labels[g] + f'_DOY{doy:03d}.nc4'
                    write_latlon_nc4(fname_out, data2d[:, :, None], [pentad], [start_time], [end_time], grid)

                if is_pentad_center:
                    if print_all_pentads:
                        data_out[g, :, :, pentad-1] = data2d
                        start_time_p[pentad-1]      = start_time
                        end_time_p[  pentad-1]      = end_time

                    if print_each_pentad:
                        fname_out = fname_base_p + '_sp_' + species_labels[g] + f'_p{pentad:02d}.nc4'
                        write_latlon_nc4(fname_out, data2d[:, :, None], [pentad], [start_time], [end_time], grid)

            for arr, empty in ((o_sum, 0.), (o_sum2, 0.), (m_sum, 0.), (m_sum2, 0.), (n_data, 0.),
                               (m_min, np.inf), (m_max, -np.inf)):
                arr[:, :, :-1] = arr[:, :, 1:]
                arr[:, :,  -1] = empty

    print(f'Read {n_read} ObsFcstAna files')
    for g in range(n_groups):
        print(f'  species {species_labels[g]}: {n_obs[g]} obs used')

    if np.any(n_obs == 0):
        raise ValueError('no obs used for species ' +
                         ', '.join(species_labels[g] for g in range(n_groups) if n_obs[g] == 0) +
                         ' (all-pentads file not written)')

    if print_all_pentads:
        for g in range(n_groups):
            fname_out = fname_base_d + '_sp_' + species_labels[g] + '_all_pentads.nc4'
            write_latlon_nc4(fname_out, data_out[g], np.arange(1, n_pentads+1), start_time_p, end_time_p, grid)

# ----------------------------------------------------------------------------
#
# statistics over the moving window (sums over the last dimension); cells with fewer
# than ndata_min obs are set to NaN

def get_window_stats(o_sum, o_sum2, m_sum, m_sum2, n_data, m_min, m_max, ndata_min):

    N = n_data.sum(axis=1)

    def mean(x):
        return np.divide(x.sum(axis=1), N, out=np.full(N.shape, np.nan), where=N > 0)

    o_mean = mean(o_sum)
    m_mean = mean(m_sum)
    o_std  = np.sqrt(np.maximum(mean(o_sum2) - o_mean**2, 0.))
    m_std  = np.sqrt(np.maximum(mean(m_sum2) - m_mean**2, 0.))

    mn = m_min.min(axis=1)
    mx = m_max.max(axis=1)
    mn[np.isinf(mn)] = np.nan
    mx[np.isinf(mx)] = np.nan

    data2d = np.array([o_mean, o_std, m_mean, m_std, N, mn, mx])

    data2d[:, N < ndata_min] = np.nan

    return data2d

# ----------------------------------------------------------------------------
#
# write statistics to nc4 file
#
# data: (n_fields, n_cells, n_pentad), with grid cell index = j*n_lon + i
#
# Variables o_mean, o_std, m_mean, m_std, n_data are (pentad, lon, lat); m_min and m_max
# are (lon, lat), the min/max over all pentads.  Each pentad is one chunk, because
# GEOSldas reads one pentad at a time.

def write_latlon_nc4(fname, data, pentads, start_time, end_time, grid):

    n_lon, n_lat, ll_lon, ll_lat, d_lon, d_lat = grid

    n_pentad = len(pentads)

    # (field, cell, pentad) -> (field, pentad, lon, lat)
    data = np.asarray(data, dtype=float).reshape(n_fields, n_lat, n_lon, n_pentad).transpose(0, 3, 2, 1)

    with Dataset(fname, 'w', format='NETCDF4') as ncid:

        ncid.createDimension('pentad', n_pentad)
        ncid.createDimension('lon',    n_lon)
        ncid.createDimension('lat',    n_lat)

        ncid.createVariable('version', 'i4')[:] = 0

        for name, value, long_name, units, axis in [
                ('ll_lon', ll_lon, 'longitude of lower left corner', 'degrees_east',  'X'),
                ('ll_lat', ll_lat, 'latitude of lower left corner',  'degrees_north', 'Y'),
                ('d_lon',  d_lon,  'longitude grid spacing',         'degrees',       'X'),
                ('d_lat',  d_lat,  'latitude grid spacing',          'degrees',       'Y')]:
            var = ncid.createVariable(name, 'f8')
            var[:] = value
            var.setncatts({'standard_name': long_name, 'long_name': long_name, 'units': units, 'axis': axis})

        var = ncid.createVariable('pentad', 'i4', ('pentad',))
        var[:] = pentads
        var.setncatts({'standard_name': 'pentad', 'long_name': 'pentad', 'units': '1', 'axis': 'T'})

        ncid.createVariable('start_time', 'f8', ('pentad',))[:] = days_since_1950(start_time)
        ncid.createVariable('end_time',   'f8', ('pentad',))[:] = days_since_1950(end_time)

        obs_units   = 'Degree of saturation (0-1)'
        model_units = 'Surface soil moisture (m^3 m^-3)'
        window      = 'for pentad calculated over all years for window length'

        fields = [
            ('o_mean', 'observation mean',               'Observation mean '               + window, obs_units  ),
            ('o_std',  'observation standard deviation', 'Observation standard deviation ' + window, obs_units  ),
            ('m_mean', 'model mean',                     'Model mean '                     + window, model_units),
            ('m_std',  'model standard deviation',       'Model standard deviation '       + window, model_units),
            ('n_data', 'number of data points',          'Number of data points '          + window, '1'        ),
            ('m_min',  'model minimum',                  'Model minimum calculated over all years',  model_units),
            ('m_max',  'model maximum',                  'Model maximum calculated over all years',  model_units)]

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)    # all-NaN cells in nanmin/nanmax
            for k, (name, standard_name, long_name, units) in enumerate(fields):
                if name == 'm_min':
                    dims, chunks, values = ('lon', 'lat'),           (n_lon, n_lat),    np.nanmin(data[k], axis=0)
                elif name == 'm_max':
                    dims, chunks, values = ('lon', 'lat'),           (n_lon, n_lat),    np.nanmax(data[k], axis=0)
                else:
                    dims, chunks, values = ('pentad', 'lon', 'lat'), (1, n_lon, n_lat), data[k]
                var = ncid.createVariable(name, 'f8', dims, chunksizes=chunks, zlib=True, complevel=5)
                var.setncatts({'standard_name': standard_name, 'long_name': long_name, 'units': units})
                var[:] = values

# ====================== EOF =========================================================
