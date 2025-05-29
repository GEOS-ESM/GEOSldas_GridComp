import os
import numpy as np

# convert data from 1-d array ("compressed") format to 2-d array (grid) format based on lat/lon coordinates
#
# 1-d input data can have one additional dimension (e.g., time)
#
# 2-d output data is lat-by-lon[-by-time]

def remap_1d_to_2d( data_1d, *, lat_1d, lon_1d):
    
    # Extract unique coordinates 
    unique_lats, indY = np.unique(lat_1d, return_inverse=True)
    unique_lons, indX = np.unique(lon_1d, return_inverse=True)

    ny = len(unique_lats)
    nx = len(unique_lons)
    
    if   data_1d.ndim == 2:                          #[n_1d, ntime]
        ntime = data_1d.shape[1]
        data_2d = np.full([ny, nx, ntime], np.nan)
        data_2d[indY, indX, :] = data_1d
    elif data_1d.ndim == 1:                          #[n_1d]
        data_2d = np.full([ny, nx],        np.nan)
        data_2d[indY, indX   ] = data_1d

    return data_2d, unique_lats, unique_lons


