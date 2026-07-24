

import os
import numpy as np

def tile2grid(tile_data, tile_coord, tile_grid):

    """
    Map (1d) tile-space data onto (2d) grid that is associated with the tile space.
    
    Parameters:
    ----------
    tile_data  : Array in tile space, shape (N_tile, N_fields)
    tile_coord : Dictionary containing tile coordinate information
    tile_grid  : Dictionary containing tile grid information

    Returns:
    -------
    grid_data : Array in grid space, shape (N_lat, N_lon, N_fields)
    """

    # Verify input datasize
    if tile_data.shape[0] != tile_coord['N_tile']:
        print(f'Error: size of tile2grid input data does not match that of N_tile')
        sys.exit()
        
    # if tile_data is 1-D [N_tile], expand into 2-D [N_tile,1]
    if tile_data.ndim == 1:
        tile_data = np.expand_dims(tile_data, axis=1)
    
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
            if  ~np.isnan(tile_data[n, k]):
                # Accumulate weighted data
                grid_data[j, i, k] += w * tile_data[n,k]
                wgrid[    j, i]    += w

        # Normalize and set no-data values
        for i in range(tile_grid['N_lon']):
            for j in range(tile_grid['N_lat']):
                if wgrid[j, i] > 0.0:
                    # Normalize by total weight
                    grid_data[j, i, k] = grid_data[j, i, k] / wgrid[j, i]
                else:
                    # Set no-data value
                    grid_data[j, i, k] = np.nan

    return grid_data.squeeze()

# ============ EOF ===============================================
