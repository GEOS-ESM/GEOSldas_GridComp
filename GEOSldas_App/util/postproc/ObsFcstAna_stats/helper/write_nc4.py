
# collection of write routines for ObsFcstAna sums and stats

from netCDF4 import Dataset, date2num

import numpy as np

# --------------------------------------------------------------------------------

def write_sums_nc4(file_path, N_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum, obs_param):
    
    nc = Dataset(file_path,'w',format='NETCDF4')
    
    tile    = nc.createDimension('tile',    N_data.shape[0])
    species = nc.createDimension('species', N_data.shape[1])

    data = nc.createVariable('obs_param_assim', 'c', ('species', ),         zlib=True, complevel=4)
    for i in range(len(obs_param)):
        data[i] = obs_param[i]['assim'] 

    data = nc.createVariable('N_data',          'i4', ('tile','species', ), zlib=True, complevel=4)
    data[:,:] = N_data

    data = nc.createVariable('obsxfcst_sum',    'f4', ('tile','species', ), zlib=True, complevel=4)
    data[:,:] = oxf_sum 

    data = nc.createVariable('obsxana_sum',     'f4', ('tile','species', ), zlib=True, complevel=4)
    data[:,:] = oxa_sum 

    data = nc.createVariable('fcstxana_sum',    'f4', ('tile','species', ), zlib=True, complevel=4)
    data[:,:] = fxa_sum 

    for key, value in data_sum.items():
        varname = key+'_sum'
        data = nc.createVariable(varname,       'f4', ('tile','species', ), zlib=True, complevel=4)
        data[:,:] = value

    for key, value in data2_sum.items():
        varname = key+'2_sum'
        data = nc.createVariable(varname,       'f4', ('tile','species', ), zlib=True, complevel=4)
        data[:,:] = value

    nc.close()

# --------------------------------------------------------------------------------

def write_stats_nc4(file_path, stats):

    nc = Dataset(file_path,'w',format='NETCDF4')
    
    tile    = nc.createDimension('tile',    stats['N_data'].shape[0])
    species = nc.createDimension('species', stats['N_data'].shape[1])

    for key, value in stats.items():
        data = nc.createVariable(key,'f4',('tile','species', ), \
                                fill_value=-9999.0, zlib=True, complevel=4)
        data[:,:] = value

    nc.close()


# ======================= EOF ========================================================= 
