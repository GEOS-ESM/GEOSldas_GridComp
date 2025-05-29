import os
import numpy as np

# convert data from an array format to a grid format based on lat/lon coordinates
def tile2grid(tile_data, *, lat_1d, lon_1d):
    
    # Extract unique coordinates 
    unique_lats, indY = np.unique(lat_1d, return_inverse=True)
    unique_lons, indX = np.unique(lon_1d, return_inverse=True)
    ny = len(unique_lats)
    nx = len(unique_lons)
    
    if   tile_data.ndim == 2:                     #[ntile, ntime]
        ntime = tile_data.shape[1]
        grid_data = np.full([ny, nx, ntime], np.nan)
        grid_data[indY, indX, :] = tile_data
    elif tile_data.ndim == 1:                     #[ntile]
        grid_data = np.full([ny, nx], np.nan)
        grid_data[indY, indX   ] = tile_data

    return grid_data, unique_lats, unique_lons


