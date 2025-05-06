from netCDF4 import Dataset, date2num
import numpy as np

def write_sums_nc4(file_path, N_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum, obs_param):
    
    nc = Dataset(file_path,'w',formate='NETCDF4')
    tile = nc.createDimension('tile', N_data.shape[0])
    species = nc.createDimension('species', N_data.shape[1])

    data = nc.createVariable('obs_param_assim', 'c', ('species', ), zlib=True, complevel=4)
    for i in range(len(obs_param)):
        data[i] = obs_param[i]['assim'] 

    data = nc.createVariable('N_data', 'i4', ('tile','species', ), zlib=True, complevel=4)
    data[:,:] = N_data

    data = nc.createVariable('obsxfcst_sum', 'f4', ('tile','species', ), zlib=True, complevel=4)
    data[:,:] = oxf_sum 

    data = nc.createVariable('obsxana_sum', 'f4', ('tile','species', ), zlib=True, complevel=4)
    data[:,:] = oxa_sum 

    data = nc.createVariable('fcstxana_sum', 'f4', ('tile','species', ), zlib=True, complevel=4)
    data[:,:] = fxa_sum 

    for key, value in data_sum.items():
        varname = key+'_sum'
        data = nc.createVariable(varname,'f4',('tile','species', ), zlib=True, complevel=4)
        data[:,:] = value

    for key, value in data2_sum.items():
        varname = key+'2_sum'
        data = nc.createVariable(varname,'f4',('tile','species', ), zlib=True, complevel=4)
        data [:,:]= value

    nc.close()

def write_stats_nc4(file_path, stats):

    nc = Dataset(file_path,'w',formate='NETCDF4')
    tile = nc.createDimension('tile', stats['N_data'].shape[0])
    species = nc.createDimension('species', stats['N_data'].shape[1])

    for key, value in stats.items():
        data = nc.createVariable(key,'f4',('tile','species', ), \
                                fill_value=-9999.0, zlib=True, complevel=4)
        data[:,:] = value

    nc.close()

def write_omf_stats_nc4(file_path, stats):
    """
    Write monthly OmF/OmA statistics to a NetCDF4 file.

    Parameters:
    -----------
    file_path : str
        Path to the output NetCDF file.

    stats : dict
        Dictionary containing time-series statistics to write. Expected keys:
            - 'time' : array-like of datetime.datetime objects (length = time)
            - All other keys (e.g., 'OmF_mean', 'OmF_stdv', etc.) should be
              3D arrays of shape (time, tile, species).

    Notes:
    ------
    - The 'time' variable is written as "days since 1900-01-01" using a standard calendar.
    - All data variables are written as float32 with zlib compression and a fill value of -9999.0.
    - The NetCDF file uses unlimited time dimension to support time series growth.
    """
        
    nc = Dataset(file_path, 'w', format='NETCDF4')
    
    # Create dimensions
    nc.createDimension('time', None)  # Unlimited dimension
    nc.createDimension('tile', stats['OmF_mean'].shape[1])
    nc.createDimension('species', stats['OmF_mean'].shape[2])
    
    # Create variables with time dimension
    for key, value in stats.items():
        if key == 'time':
            var = nc.createVariable(key, 'f8', ('time',), zlib=True)
            var.units = 'days since 1900-01-01'
            var.calendar = 'standard'
            var[:] = date2num(value, units=var.units, calendar=var.calendar)
        else:
            var = nc.createVariable(key, 'f4', ('time', 'tile', 'species'), 
                                  fill_value=-9999.0, zlib=True, complevel=4)
            var[:] = value
    
    nc.close()    


def write_omf_grouped_stats_nc4(file_path, group_stats):
    """
    Write sensor-grouped, per-tile OMF/OMA statistics to NetCDF.
    
    Parameters:
        file_path   : str
        group_stats : dict
            {
              'SMAP': {'OmF_mean': (time, tile), ...},
              'SMOS': {...},
              ...
              'time': [datetime1, datetime2, ...]
            }
    """
    # Prepare metadata
    group_names = [g for g in group_stats if g != 'time']
    n_group = len(group_names)
    n_time = len(group_stats['time'])
    n_tile = next(iter(group_stats[group_names[0]]['OmF_mean'])).shape[0] if n_time == 1 else group_stats[group_names[0]]['OmF_mean'].shape[1]

    # Create file
    nc = Dataset(file_path, 'w', format='NETCDF4')

    # Create dimensions
    nc.createDimension('time', n_time)
    nc.createDimension('tile', n_tile)
    nc.createDimension('group', n_group)

    # Write time
    time_var = nc.createVariable('time', 'f8', ('time',))
    time_var.units = 'days since 1900-01-01'
    time_var.calendar = 'standard'
    time_var[:] = date2num(group_stats['time'], units=time_var.units, calendar=time_var.calendar)

    # Write group names as string variable
    group_name_len = max(len(name) for name in group_names)
    nc.createDimension('name_strlen', group_name_len)
    group_name_var = nc.createVariable('group_name', 'S1', ('group', 'name_strlen'))
    for i, name in enumerate(group_names):
        group_name_var[i, :len(name)] = np.array(list(name), 'S1')

    # Variables to write
    metrics = ['N_data', 'OmF_mean', 'OmF_stdv', 'OmF_norm_mean', 'OmF_norm_stdv', 'OmA_mean', 'OmA_stdv']
    for metric in metrics:
        var = nc.createVariable(metric, 'f4', ('time', 'tile', 'group'), fill_value=-9999.0, zlib=True)

        # Build (time, tile, group) array
        data = np.full((n_time, n_tile, n_group), np.nan, dtype=np.float32)
        for g_idx, g_name in enumerate(group_names):
            data[:, :, g_idx] = group_stats[g_name][metric]
        var[:] = data

    nc.close()    
