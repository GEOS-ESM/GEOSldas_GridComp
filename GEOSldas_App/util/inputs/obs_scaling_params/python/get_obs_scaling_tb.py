
# Compute climatological mean and stdv of Tb observations and model forecasts from
# GEOSldas ObsFcstAna files in tile space, for use in the z-score scaling of Tb
# observations before assimilation (GEOSldas scale_obs_Tb_zscore).
#
# For each day of the climatological year, the statistics are computed from all
# obs/fcst pairs within a moving window of w_days days (across all years), and
# stored at the center of the window.  One file is written for each pentad (and,
# optionally, for each DOY) in the big-endian Fortran sequential binary format
# read by GEOSldas.  One orbit (A or D) per call; H and V pol for each angle.
#
# Python version of get_model_and_obs_clim_stats.m (for hscale=0).
#
# amfox - Oct 2026

import os
import sys;        sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../shared/python'))

import numpy as np

from datetime          import datetime, timedelta

from read_GEOSldas      import read_ObsFcstAna_nc4, read_tilecoord, read_tilegrids
from EASEv2            import EASEv2_latlon2ind, EASEv2_ind2latlon
from obs_scaling_utils import get_ObsFcstAna_fname, get_doy_pentad_range, days_in_month, nodata
from obs_scaling_utils import any_ObsFcstAna_file

n_fields = 14

def get_obs_scaling_tb(exp_path, exp_run, domain, obs_param, descr_contains, orbit, angles,
                       run_months, first_year, last_year, dt_assim, t0_assim,
                       w_days, ndata_min, prefix, out_dir, convert_grid=True, print_each_doy=False):

    # orbit       : 1 (ascending) or 2 (descending)
    # convert_grid: True  = use EASEv2 M36 "administering" tiles (one tile per M36 grid cell)
    #               False = use all model tiles
    #
    # The other arguments are checked in run_obs_scaling_params.py.

    if orbit not in (1, 2):
        raise ValueError('orbit must be 1 (ascending) or 2 (descending)')

    inpath  = exp_path + '/' + exp_run + '/output/' + domain
    outpath = inpath + '/stats/' + out_dir

    t0_assim = t0_assim % dt_assim
    cycles   = range(t0_assim, 86400, dt_assim)

    if not any_ObsFcstAna_file(inpath, exp_run, run_months, first_year, last_year, cycles):
        raise ValueError(f'no ObsFcstAna files for the requested period in {inpath}/ana/ens_avg')

    os.makedirs(outpath, exist_ok=True)

    tile_coord   = read_tilecoord(inpath + '/rc_out/' + exp_run + '.ldas_tilecoord.bin')
    _, tile_grid = read_tilegrids(inpath + '/rc_out/' + exp_run + '.ldas_tilegrids.bin')

    N_tile = tile_coord['N_tile']

    # output tiles, and index into output tiles for each model tile (-1 if not an output tile)

    if convert_grid:
        out_tiles = get_administering_tiles(tile_coord, tile_grid)
    else:
        out_tiles = np.arange(N_tile)
    N_out = len(out_tiles)

    model_to_out            = np.full(N_tile, -1)
    model_to_out[out_tiles] = np.arange(N_out)

    print(f"Selected {N_out:,} output tiles from {N_tile:,} model tiles")

    # Tb species: {species: (pol index, angle index)}

    tb_species = select_tb_species(obs_param, descr_contains, orbit, angles)

    # obs lat/lon must be close to the lat/lon of their tile (M36 obs on M36 tiles: exact)

    is_m36_native = 'M36' in domain and all('_E' not in descr for descr in tb_species_descr(obs_param, tb_species))
    tol = 1e-3 if is_m36_native else 2.

    # sums for each pol, output tile, angle and day of the moving window:
    #  sum(obs), sum(obs^2), sum(fcst), sum(fcst^2), N

    ring = np.zeros((5, 2, N_out, len(angles), w_days))

    # output file names

    d1, d2, p1, p2 = get_doy_pentad_range(run_months, first_year, last_year)

    orb_tag      = 'A' if orbit == 1 else 'D'
    fname_base_d = (outpath + '/' + prefix + f'{min(first_year)}_doy{d1}_{max(last_year)}_doy{d2}' +
                    f'_hscale_0.00_W_{w_days}d_Nmin_{ndata_min}_{orb_tag}')
    fname_base_p = (outpath + '/' + prefix + f'{min(first_year)}_p{p1}_{max(last_year)}_p{p2}' +
                    f'_hscale_0.00_W_{round(w_days / 5)}p_Nmin_{ndata_min}_{orb_tag}')

    n_read = 0                          # number of ObsFcstAna files read
    n_obs  = 0                          # number of Tb obs used

    count = 0

    for imonth, month in enumerate(run_months):
        for day in range(1, days_in_month[month-1] + 1):

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
                    valid = np.isfinite(ofa['obs_fcst']) & np.isfinite(ofa['obs_obs'])

                    for species, (ipol, iang) in tb_species.items():

                        mask = valid & (ofa['obs_species'] == species)
                        if not np.any(mask):
                            continue

                        tiles = ofa['obs_tilenum'][mask].astype(np.int64) - 1

                        if np.any(tiles < 0) or np.any(tiles >= N_tile):
                            raise ValueError(f'species {species}: tile number out of range in {fname}')
                        if np.unique(tiles).size != tiles.size:
                            raise ValueError(f'species {species}: more than one obs per tile in {fname}')

                        out = model_to_out[tiles]
                        if np.any(out < 0):
                            raise ValueError(f'species {species}: obs on {np.count_nonzero(out < 0)} '
                                             f'tiles that are not output tiles in {fname}')

                        if (np.any(np.abs(tile_coord['com_lon'][tiles] - ofa['obs_lon'][mask]) > tol) or
                            np.any(np.abs(tile_coord['com_lat'][tiles] - ofa['obs_lat'][mask]) > tol)):
                            raise ValueError(f'obs lat/lon and tile lat/lon do not match in {fname}')

                        obs  = ofa['obs_obs' ][mask].astype(np.float64)
                        fcst = ofa['obs_fcst'][mask].astype(np.float64)
                        n_obs += len(out)

                        np.add.at(ring[0, ipol, :, iang, idx], out, obs    )
                        np.add.at(ring[1, ipol, :, iang, idx], out, obs**2 )
                        np.add.at(ring[2, ipol, :, iang, idx], out, fcst   )
                        np.add.at(ring[3, ipol, :, iang, idx], out, fcst**2)
                        np.add.at(ring[4, ipol, :, iang, idx], out, 1.     )

            if count < w_days:
                continue

            # The window is full: compute the statistics, write them for the center of
            # the window, and shift the window by one day.

            end_time   = datetime(2014, month, day, hour, minute)
            start_time = end_time - timedelta(seconds=w_days*86400)
            doy        = (end_time - timedelta(seconds=w_days*86400//2)).timetuple().tm_yday

            is_pentad_center = (doy + 2) % 5 == 0

            fnames_out = []
            if print_each_doy:
                fnames_out.append(fname_base_d + f'_DOY{doy:03d}.bin')
            if is_pentad_center:
                fnames_out.append(fname_base_p + f'_p{(doy + 2) // 5:02d}.bin')

            if fnames_out:

                data_out = get_window_stats_tb(ring, ndata_min, w_days - w_days//2 - 1)

                # 2nd header entry is 0 (as written by the MATLAB scripts), ndata_min is in the file name
                for fname_out in fnames_out:
                    write_tb_scaling_bin(fname_out, 1 if orbit == 1 else 0, 0, start_time, end_time, angles,
                                         tile_coord['com_lon'][out_tiles], tile_coord['com_lat'][out_tiles],
                                         tile_coord['tile_id'][out_tiles], data_out)

            ring[..., :-1] = ring[..., 1:]
            ring[...,  -1] = 0.

    print(f'Read {n_read} ObsFcstAna files, {n_obs} Tb obs used')

    if n_obs == 0:
        raise ValueError(f'no Tb obs used; the files written to {outpath} contain no data')

# ----------------------------------------------------------------------------
#
# Tb species for orbit and angles: {species: (pol index, angle index)}
# (exactly one H-pol and one V-pol species for each angle)

def select_tb_species(obs_param, descr_contains, orbit, angles):

    tb_species = {}

    for pol in (1, 2):
        for iang, angle in enumerate(angles):
            matches = [p for p in obs_param
                       if p['varname'] == 'Tb' and descr_contains in p['descr'] and
                       p['orbit'] == orbit and p['pol'] == pol and
                       np.any(np.isclose(p['ang'], angle, atol=1e-6, rtol=0.))]
            if len(matches) != 1:
                raise ValueError(f'expected one Tb species for orbit={orbit}, pol={pol}, angle={angle}; '
                                 f'found {[p["descr"] for p in matches]}')
            tb_species[int(matches[0]['species'])] = (pol - 1, iang)
            print(f"Tb species {int(matches[0]['species'])}: {matches[0]['descr']}")

    return tb_species

def tb_species_descr(obs_param, tb_species):

    return [p['descr'] for p in obs_param if int(p['species']) in tb_species]

# ----------------------------------------------------------------------------
#
# statistics over the moving window, (14, N_out, N_ang):
#
#   0- 4: H-pol o_mean, o_std, m_mean, m_std, N   (nodata if N < ndata_min, or, for a
#                                                  95-day window, if the 31 center days
#                                                  have fewer than ndata_min/2.5 obs)
#   5- 9: V-pol o_mean, o_std, m_mean, m_std, N
#  10-13: H-pol o_mean, m_mean, V-pol o_mean, m_mean for the center day of the window

def get_window_stats_tb(ring, ndata_min, i_center):

    data_out = np.full((n_fields, ring.shape[2], ring.shape[3]), nodata, dtype=np.float32)

    for ipol in range(2):

        N = ring[4, ipol].sum(axis=-1)

        def mean(x):
            return np.divide(x, N, out=np.full(N.shape, np.nan), where=N > 0)

        o_mean = mean(ring[0, ipol].sum(axis=-1))
        m_mean = mean(ring[2, ipol].sum(axis=-1))
        o_std  = np.sqrt(np.maximum(mean(ring[1, ipol].sum(axis=-1)) - o_mean**2, 0.))
        m_std  = np.sqrt(np.maximum(mean(ring[3, ipol].sum(axis=-1)) - m_mean**2, 0.))

        valid = N >= ndata_min
        if ring.shape[-1] == 95:
            N_inner = ring[4, ipol, :, :, i_center-15:i_center+16].sum(axis=-1)
            valid  &= N_inner >= ndata_min / 2.5
        for k, value in enumerate((o_mean, o_std, m_mean, m_std, N)):
            data_out[5*ipol + k] = np.where(valid, value, nodata)

        N_c = ring[4, ipol, :, :, i_center]
        o_c = np.divide(ring[0, ipol, :, :, i_center], N_c, out=np.full(N_c.shape, np.nan), where=N_c > 0)
        m_c = np.divide(ring[2, ipol, :, :, i_center], N_c, out=np.full(N_c.shape, np.nan), where=N_c > 0)

        data_out[10 + 2*ipol] = np.where(N_c > 0, o_c, nodata)
        data_out[11 + 2*ipol] = np.where(N_c > 0, m_c, nodata)

    return data_out

# ----------------------------------------------------------------------------
#
# write Tb scaling parameters (big-endian Fortran sequential binary, as read by
# GEOSldas scale_obs_Tb_zscore and read_Tb_scale_param.py)

def write_tb_scaling_bin(fname, asc_flag, header2, start_time, end_time, angles, lon, lat, tile_id, data):

    def write_record(fid, values, dtype):
        payload = np.asarray(values, dtype=dtype).tobytes()
        tag     = np.array([len(payload)], dtype='>i4').tobytes()
        fid.write(tag + payload + tag)

    def time_values(t):
        return [t.year, t.month, t.day, t.hour, t.minute]

    N_out = len(lon)
    N_ang = len(angles)

    with open(fname, 'wb') as fid:
        write_record(fid, [asc_flag, header2],   '>i4')
        write_record(fid, time_values(start_time), '>i4')
        write_record(fid, time_values(end_time),   '>i4')
        write_record(fid, [N_out, N_ang, 1],      '>i4')
        write_record(fid, angles,                 '>f4')
        write_record(fid, lon,                    '>f4')
        write_record(fid, lat,                    '>f4')
        write_record(fid, tile_id,                '>i4')
        for k in range(n_fields):
            for iang in range(N_ang):
                write_record(fid, data[k, :, iang], '>f4')

# ----------------------------------------------------------------------------
#
# EASEv2 M36 "administering" tiles: for each M36 grid cell that contains a model tile,
# the model tile at the (slightly shifted) center of the M36 cell, unless the cell
# center is outside the bounding box of that tile AND more than fov_km from its center.
# Returns indices (0-based) into the model tiles.

def get_administering_tiles(tile_coord, tile_grid, shift_lon=0.01, shift_lat=0.005, fov_km=20.):

    if tile_grid['gridtype'].replace('-', '_').endswith('M36'):
        model_grid = 'M36'
    elif tile_grid['gridtype'].replace('-', '_').endswith('M09'):
        model_grid = 'M09'
    else:
        raise ValueError(f"tile grid {tile_grid['gridtype']} not supported")

    # M36 grid cells that contain a model tile, and their (shifted) centers

    m36_row, m36_col = EASEv2_latlon2ind(tile_coord['com_lat'], tile_coord['com_lon'], 'M36')
    m36_cells        = np.unique(np.column_stack((m36_row, m36_col)), axis=0)
    lat, lon         = EASEv2_ind2latlon(m36_cells[:, 0], m36_cells[:, 1], 'M36')
    lat              = lat + shift_lat
    lon              = lon + shift_lon

    # model tile for each cell of the (domain) model tile grid

    i_tile = tile_coord['i_indg'].astype(np.int64) - tile_grid['i_offg'] - tile_grid['ind_base']
    j_tile = tile_coord['j_indg'].astype(np.int64) - tile_grid['j_offg'] - tile_grid['ind_base']
    if (np.any(i_tile < 0) or np.any(i_tile >= tile_grid['N_lon']) or
        np.any(j_tile < 0) or np.any(j_tile >= tile_grid['N_lat'])):
        raise ValueError('tilecoord i_indg/j_indg outside the domain tile grid')
    cell_to_tile = {int(c): n for n, c in enumerate(j_tile * tile_grid['N_lon'] + i_tile)}
    if len(cell_to_tile) != len(i_tile):
        raise ValueError('more than one tile in a grid cell')

    # model tile at each (shifted) M36 cell center

    row, col = EASEv2_latlon2ind(lat, lon, model_grid)
    i_cell   = col - tile_grid['i_offg'] - tile_grid['ind_base']
    j_cell   = row - tile_grid['j_offg'] - tile_grid['ind_base']

    dist_y = fov_km * 180. / np.pi / 6371.
    dist_x = dist_y / np.cos(np.deg2rad(lat))

    tiles = []
    for n in range(len(lat)):
        if i_cell[n] < 0 or i_cell[n] >= tile_grid['N_lon'] or j_cell[n] < 0 or j_cell[n] >= tile_grid['N_lat']:
            continue
        t = cell_to_tile.get(int(j_cell[n] * tile_grid['N_lon'] + i_cell[n]))
        if t is None:
            continue
        outside_bbox = (lon[n] < tile_coord['min_lon'][t] or lon[n] > tile_coord['max_lon'][t] or
                        lat[n] < tile_coord['min_lat'][t] or lat[n] > tile_coord['max_lat'][t])
        too_far      = (abs(lon[n] - tile_coord['com_lon'][t]) > dist_x[n] or
                        abs(lat[n] - tile_coord['com_lat'][t]) > dist_y)
        if not (outside_bbox and too_far):
            tiles.append(t)

    tiles = np.array(tiles, dtype=np.int64)
    if np.unique(tiles).size != tiles.size:
        raise ValueError('administering-tile selection gave duplicate tiles')

    return tiles

# ====================== EOF =========================================================
