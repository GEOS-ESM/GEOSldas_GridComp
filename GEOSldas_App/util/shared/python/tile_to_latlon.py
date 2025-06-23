import os
import numpy as np

# convert data from 1-d array ("compressed") format to 2-d array (grid) format based on lat/lon coordinates
#
# 1-d input data can have one additional dimension (e.g., time)
#
# 2-d output data is lat-by-lon[-by-time]

def get_grid_resolution(fov_radius=0.5):
    '''
    Estimate the resolution in degrees for the target lat/lon grid when remapping 1-D
    tile data to regular 2-D lat/lon grid for plotting purposes.
    The resolution is an approximation based on the observation FOV radius information. 
    
    '''    
    return np.round(fov_radius * 2., decimals=3)
        
def interpolate_to_latlon(data_1d, lat, lon, resolution):
    """
    For cubedsphere grid: use 
    
    """
    from scipy.interpolate import griddata
    from scipy.spatial.distance import cdist
    
    lat_grid = np.arange(-90+resolution/2, 90, resolution)
    lon_grid = np.arange(-180+resolution/2, 180, resolution)
    lon_2d, lat_2d = np.meshgrid(lon_grid, lat_grid)

    ny = len(lat_grid)
    nx = len(lon_grid)

    nf = data_1d.shape[-1]
    data_2d = np.full([ny, nx, nf], np.nan)
    for i in range(nf):
        data = data_1d[:,i]
        data_2d[:,:,i] = griddata((lat.flatten(), lon.flatten()), data.flatten(), (lat_2d, lon_2d), method='linear')
     
    return data_2d, lat_2d, lon_2d

def remap_1d_to_2d( data_1d,  lat_1d, lon_1d):

    """
    For EASE grid: use remapping
    
    """
    
    # Extract unique coordinates 
    unique_lats, indY = np.unique(lat_1d, return_inverse=True)
    unique_lons, indX = np.unique(lon_1d, return_inverse=True)

    ny = len(unique_lats)
    nx = len(unique_lons)
    
    nf = data_1d.shape[1]
    data_2d = np.full([ny, nx, nf], np.nan)
    data_2d[indY, indX, :] = data_1d
        
    lon_2d,lat_2d = np.meshgrid(unique_lons, unique_lats)
    
    return data_2d, lat_2d, lon_2d

def tile2grid(tile_data, tile_coord, tile_grid):

    N_fields = tile_data.shape[-1]
    # Initialize grid data array
    grid_data = np.zeros((tile_grid['N_lat'], tile_grid['N_lon'], N_fields))

    for k in range(N_fields):
        # Initialize weight grid for current field
        wgrid = np.zeros((tile_grid['N_lat'], tile_grid['N_lon']))
        
        # Loop through tile space
        for n in range(tile_coord['N_tile']):
            # Calculate grid indices (adjust for Python's 0-based indexing)
            i = tile_coord['i_indg'][n] - (tile_grid['i_offg'] - (1-tile_grid['ind_base'])) - 1
            j = tile_coord['j_indg'][n] - (tile_grid['j_offg'] - (1-tile_grid['ind_base'])) - 1
            
            # Get weight
            w = tile_coord['frac_cell'][n]
            
            # Check if current tile data is valid (not no-data)
            if  ~np.isnan(tile_data[k, n]):
                # Accumulate weighted data
                grid_data[j, i, k] += w * tile_data[n,k]
                wgrid[j, i] += w
        
        # Normalize and set no-data values
        for i in range(tile_grid['N_lon']):
            for j in range(tile_grid['N_lat']):
                if wgrid[j, i] > 0.0:
                    # Normalize by total weight
                    grid_data[j, i, k] = grid_data[j, i, k] / wgrid[j, i]
                else:
                    # Set no-data value
                    grid_data[j, i, k] = np.nan
    
    lat_2d = np.arange(tile_grid['ll_lat']+0.5*tile_grid['dlat'], tile_grid['ur_lat'], tile_grid['dlat'])
    lon_2d = np.arange(tile_grid['ll_lon']+0.5*tile_grid['dlon'], tile_grid['ur_lon'], tile_grid['dlon'])
    
def tile_to_latlon( tile_data, tile_coord, tile_grid, resolution, nodata=1.e15, nodata_tol=1.e-4):
    
    """
    Convert tile-space data to grid-space data.
    
    Parameters:
    ----------
    tile_data : Array in tile space, shape (N_tile, N_fields)
    tile_coord : Dictionary containing tile coordinate information
    tile_grid : Dictionary containing grid information
    resolution : target grid resolution
    nodata : Value for missing data
    nodata_tol : Tolerance for comparing values to nodata
        
    Returns:
    -------
    grid_data : Array in regular grid space, shape (N_lat, N_lon, N_fields)
    lat_2d: Lat array in regular grid space, shape (N_lat, N_lon)
    lon_2d: Lon array in regular grid space, shape (N_lat, N_lon)
    
    """
    tile_data[np.abs(tile_data - nodata) < nodata_tol]  = np.nan
    
    # Verifity input datasize
    if tile_data.shape[0] != tile_coord['N_tile']:
        print(f'Error: tile2grid input data size don''t match N_tile')
        sys.exit()
        
    # if tile_data is 1-D [N_tile], expand into 2-D [N_tile,1]
    if tile_data.ndim == 1:
        tile_data = np.expand_dims(tile_data, axis=1)

    nf = tile_data.shape[-1]

    lat_grid = np.arange(-90+resolution/2, 90, resolution)
    lon_grid = np.arange(-180+resolution/2, 180, resolution)
    lon_2d, lat_2d = np.meshgrid(lon_grid, lat_grid)

    grid_data = np.full([len(lat_grid), len(lon_grid), nf], np.nan)

    for v in range(nf):        
        for k in range(tile_coord['N_tile']):
            if np.any([~np.isnan(tile_data[k,:])], axis=1):
                j = np.searchsorted(lat_grid + resolution/2., tile_coord['com_lat'][k])
                i = np.searchsorted(lon_grid + resolution/2., tile_coord['com_lon'][k])
                # sometime the missing points are filled with zeros. taking largest value within
                # the target grid cell avoild picking up unwanted zeros.
                grid_data[j,i,v] = np.nanmax([tile_data[k,v], grid_data[j,i,v]])
                                            
    if grid_data.shape[-1] == 1:
        grid_data = grid_data[:,:,0]
        
    return grid_data, lat_2d, lon_2d

