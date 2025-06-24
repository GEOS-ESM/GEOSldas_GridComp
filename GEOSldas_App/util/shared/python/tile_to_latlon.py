import os
import numpy as np

# 
        
def tile_to_latlongrid( tile_data, tile_coord, resolution, nodata=1.e15, nodata_tol=1.e-4):
    
    """
    Map (1d) tile-space data onto (2d) regular lat/lon grid of arbitrary "resolution".
    
    Parameters:
    ----------
    tile_data  : Array in tile space, shape (N_tile, N_fields)
    tile_coord : Dictionary containing tile coordinate information
    resolution : target grid resolution
    nodata     : Value for missing data
    nodata_tol : Tolerance for comparing values to nodata
        
    Returns:
    -------
    grid_data : Array in regular grid space, shape (N_lat, N_lon, N_fields)
    lat_2d    : Lat array in regular grid space, shape (N_lat, N_lon)
    lon_2d    : Lon array in regular grid space, shape (N_lat, N_lon)
    
    """
    
    tile_data[np.abs(tile_data - nodata) < nodata_tol]  = np.nan
    
    # Verify input datasize
    if tile_data.shape[0] != tile_coord['N_tile']:
        print(f'Error: size of tile2grid input data does not match that of N_tile')
        sys.exit()
        
    # if tile_data is 1-D [N_tile], expand into 2-D [N_tile,1]
    if tile_data.ndim == 1:
        tile_data = np.expand_dims(tile_data, axis=1)

    nf = tile_data.shape[-1]

    lat_grid = np.arange( -90+resolution/2,  90, resolution)
    lon_grid = np.arange(-180+resolution/2, 180, resolution)

    lon_2d, lat_2d = np.meshgrid(lon_grid, lat_grid)

    grid_data = np.full([len(lat_grid), len(lon_grid), nf], np.nan)

    for v in range(nf):        
        for k in range(tile_coord['N_tile']):
            if np.any([~np.isnan(tile_data[k,:])], axis=1):
                j = np.searchsorted(lat_grid + resolution/2., tile_coord['com_lat'][k])
                i = np.searchsorted(lon_grid + resolution/2., tile_coord['com_lon'][k])
                # sometime the missing points are filled with zeros. taking largest value within
                # the target grid cell avoids picking up unwanted zeros.
                grid_data[j,i,v] = np.nanmax([tile_data[k,v], grid_data[j,i,v]])
                                            
    if grid_data.shape[-1] == 1:
        grid_data = grid_data[:,:,0]
        
    return grid_data, lat_2d, lon_2d

# =========================== EOF ============================================================
