
# collection of readers for GEOSldas output files
#
#   - read_obs_param()
#   - read_tilecoord()
#   - read_tilegrids()
#   - read_ObsFcstAna()
#   - read_ObsFcstAna_nc4()
#   - read_catparam()

import struct
import os
import numpy as np

# ----------------------------------------------------------------------------
#
# reader for GEOSldas obsparam file (ASCII)

def read_obs_param(fname):
    print(f"Reading {fname}")

    with open(fname, 'r') as fid:
        N_obs_param = int(fid.readline().strip())

        obs_param = []
        for _ in range(N_obs_param):
            param = {}
            param['descr']          =       fid.readline().strip().strip('"')
            param['species']        = float(fid.readline().strip())
            param['orbit']          = float(fid.readline().strip())
            param['pol']            = float(fid.readline().strip())

            param['N_ang']          = int(float(fid.readline().strip()))

            param['ang']            = np.array([float(x) for x in fid.readline().split()])

            param['freq']           = float(fid.readline().strip())
            param['FOV']            = float(fid.readline().strip())
            param['FOV_units']      =       fid.readline().strip().strip('"')
            param['assim']          =       fid.readline().strip()
            param['scale']          =       fid.readline().strip()
            param['getinnov']       =       fid.readline().strip()
            param['RTM_ID']         = float(fid.readline().strip())
            param['bias_Npar']      = float(fid.readline().strip())
            param['bias_trel']      = float(fid.readline().strip())
            param['bias_tcut']      = float(fid.readline().strip())
            param['nodata']         = float(fid.readline().strip())
            param['varname']        =       fid.readline().strip().strip('"')
            param['units']          =       fid.readline().strip().strip('"')
            param['fcstvarname']    =       fid.readline().strip().strip('"')
            param['fcstunits']      =       fid.readline().strip().strip('"')
            param['path']           =       fid.readline().strip().strip('"')
            param['name']           =       fid.readline().strip().strip('"')
            param['maskpath']       =       fid.readline().strip().strip('"')
            param['maskname']       =       fid.readline().strip().strip('"')
            param['scalepath']      =       fid.readline().strip().strip('"')
            param['scalename']      =       fid.readline().strip().strip('"')
            param['flistpath']      =       fid.readline().strip().strip('"')
            param['flistname']      =       fid.readline().strip().strip('"')
            param['errstd']         = float(fid.readline().strip())
            param['std_normal_max'] = float(fid.readline().strip())
            param['zeromean']       =       fid.readline().strip()
            param['coarsen_pert']   =       fid.readline().strip()
            param['xcorr']          = float(fid.readline().strip())
            param['ycorr']          = float(fid.readline().strip())
            param['adapt']          = float(fid.readline().strip())

            obs_param.append(param)

    print(f"Done reading obs_param for {N_obs_param} species")

    return obs_param

# ----------------------------------------------------------------------------
#
# reader for GEOSldas tilecoord file (binary)

def read_tilecoord(fname):
    int_precision   = 'i'
    float_precision = 'f'

    # SPECIFY endianness
    machfmt = '<'                    # '>' for big-endian, '<' for little-endian

    print(f"reading from {fname}")

    tile_coord = {}

    with open(fname, 'rb') as ifp:
        fortran_tag          = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
        tile_coord['N_tile'] = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
        fortran_tag          = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

        Nt = tile_coord['N_tile']

        fields = ['tile_id', 'typ', 'pfaf', 'com_lon', 'com_lat', 'min_lon', 'max_lon',
                  'min_lat', 'max_lat', 'i_indg', 'j_indg', 'frac_cell', 'frac_pfaf',
                  'area', 'elev']

        for field in fields:
            this_dtype        = int_precision if field in ['tile_id', 'typ', 'pfaf', 'i_indg', 'j_indg'] else float_precision
            
            fortran_tag       = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tile_coord[field] = np.frombuffer(ifp.read(Nt * 4), dtype=f'{machfmt}{this_dtype}')
            fortran_tag       = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

    print("done reading file")
    return tile_coord

# ----------------------------------------------------------------------------
#
# reader for GEOSldas tilegrids file (binary)

def read_tilegrids(fname):
    """
    Read tile grid information from file and return global and domain grid structures.
    
    Parameters:
    ----------
    fname : str
        Path to the input file (either .txt or .bin)
        
    Returns:
    -------
    tile_grid_g : dict
        Global tile grid information
    tile_grid_d : dict
        Domain tile grid information
    """
    
    # Set endian format
    machfmt = '<'            # '>' for big-endian, '<' for little-endian
    
    # Read binary file
    print(f'reading from {fname}')
    
    with open(fname, 'rb') as ifp:
        # Read "global" and "domain" records
        for grid in ['global','domain']:
            
            tile_grid = {}
            
            fortran_tag           = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tile_grid['gridtype'] = ifp.read(40).decode('ascii').strip('\x00')
            tile_grid['ind_base'] = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tile_grid['i_dir']    = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tile_grid['j_dir']    = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tile_grid['N_lon']    = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tile_grid['N_lat']    = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tile_grid['i_offg']   = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tile_grid['j_offg']   = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tile_grid['ll_lon']   = struct.unpack(f'{machfmt}f', ifp.read(4))[0]
            tile_grid['ll_lat']   = struct.unpack(f'{machfmt}f', ifp.read(4))[0]
            tile_grid['ur_lon']   = struct.unpack(f'{machfmt}f', ifp.read(4))[0]
            tile_grid['ur_lat']   = struct.unpack(f'{machfmt}f', ifp.read(4))[0]
            tile_grid['dlon']     = struct.unpack(f'{machfmt}f', ifp.read(4))[0]
            tile_grid['dlat']     = struct.unpack(f'{machfmt}f', ifp.read(4))[0]
            fortran_tag           = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            if 'global' in grid:
                tile_grid_g = tile_grid
            else:
                tile_grid_d = tile_grid
            
    return tile_grid_g, tile_grid_d

# ----------------------------------------------------------------------------
#
# reader for GEOSldas ObsFcstAna file (binary)

def read_ObsFcstAna(fname, isLDASsa=False):

    # Initialize outputs
    nodata = -9999

    date_time = {
        'year'  : nodata,
        'month' : nodata,
        'day'   : nodata,
        'hour'  : nodata,
        'min'   : nodata,
        'sec'   : nodata,
        'dofyr' : nodata,
        'pentad': nodata
    }

    obs_assim   = []
    obs_species = []
    obs_tilenum = []
    obs_lon     = []
    obs_lat     = []
    obs_obs     = []
    obs_obsvar  = []
    obs_fcst    = []
    obs_fcstvar = []
    obs_ana     = []
    obs_anavar  = []

    # SPECIFY endianness
    machfmt = '<'                    # '>' for big-endian, '<' for little-endian
    
    if os.path.exists(fname):
        print(f"reading from {fname}")

        with open(fname, 'rb') as ifp:
            # Read N_obs and time stamp entry
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            N_obs       = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            year, month, day, hour, minute, second, dofyr, pentad = struct.unpack(f'{machfmt}8i', ifp.read(32))
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            date_time = {
                'year'  : year,
                'month' : month,
                'day'   : day,
                'hour'  : hour,
                'min'   : minute,
                'sec'   : second,
                'dofyr' : dofyr,
                'pentad': pentad
            }

            # Read observation assim flag
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tmp_data = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}i').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            obs_assim = np.zeros(N_obs, dtype=int)
            obs_assim[tmp_data != 0] = 1

            # Read species information
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_species = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}i').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read tile number information
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_tilenum = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}i').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read longitude
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_lon     = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read latitude
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_lat     = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation value
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_obs     = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation variance
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_obsvar  = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space model forecast value
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_fcst    = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space model forecast variance
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_fcstvar = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space analysis value
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_ana     = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space analysis variance
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_anavar  = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # No-data check
            obs_obsvar[ obs_obsvar  == nodata] = np.nan
            obs_fcst[   obs_fcst    == nodata] = np.nan
            obs_fcstvar[obs_fcstvar == nodata] = np.nan
            obs_ana[    obs_ana     == nodata] = np.nan
            obs_anavar[ obs_anavar  == nodata] = np.nan

    else:
        print(f"file does not exist: {fname}")

    return {'date_time'  : date_time, 
            'obs_assim'  : obs_assim, 
            'obs_species': obs_species, 
            'obs_tilenum': obs_tilenum, 
            'obs_lon'    : obs_lon, 
            'obs_lat'    : obs_lat,
            'obs_obs'    : obs_obs, 
            'obs_obsvar' : obs_obsvar, 
            'obs_fcst'   : obs_fcst, 
            'obs_fcstvar': obs_fcstvar, 
            'obs_ana'    : obs_ana, 
            'obs_anavar' : obs_anavar}

# ----------------------------------------------------------------------------
#
# reader for GEOSldas ObsFcstAna file (NetCDF-4)

def read_ObsFcstAna_nc4(fname):

    from datetime import datetime
    from netCDF4 import Dataset

    nodata = -9999

    date_time = {
        'year'  : nodata,
        'month' : nodata,
        'day'   : nodata,
        'hour'  : nodata,
        'min'   : nodata,
        'sec'   : nodata,
        'dofyr' : nodata,
        'pentad': nodata
    }

    obs_assim   = []
    obs_species = []
    obs_tilenum = []
    obs_lon     = []
    obs_lat     = []
    obs_obs     = []
    obs_obsvar  = []
    obs_fcst    = []
    obs_fcstvar = []
    obs_ana     = []
    obs_anavar  = []

    if os.path.exists(fname):
        print(f"reading from {fname}")

        date_string = os.path.basename(fname).split('.')[-2].rstrip('z')
        dt = datetime.strptime(date_string, '%Y%m%d_%H%M')
        dofyr = int(dt.strftime('%j'))
        is_leap_year = (dt.year % 4 == 0 and (dt.year % 100 != 0 or dt.year % 400 == 0))
        pentad = (dofyr-1)//5 + 1
        if is_leap_year and dofyr >= 59:
            pentad = (dofyr-2)//5 + 1

        date_time = {
            'year'  : dt.year,
            'month' : dt.month,
            'day'   : dt.day,
            'hour'  : dt.hour,
            'min'   : dt.minute,
            'sec'   : dt.second,
            'dofyr' : dofyr,
            'pentad': pentad
        }

        with Dataset(fname, 'r') as ncid:
            def get_int(name):
                return np.array(ncid.variables[name][:], dtype=int)

            def get_float(name, mask_nodata=True):
                var = ncid.variables[name]
                data = var[:]
                if np.ma.isMaskedArray(data):
                    if mask_nodata:
                        data = data.filled(np.nan)
                    else:
                        data = data.filled(getattr(var, '_FillValue', np.nan))
                data = np.array(data, dtype=np.float32)

                if mask_nodata:
                    fill_value = getattr(var, '_FillValue', None)
                    if fill_value is not None:
                        data[data == fill_value] = np.nan

                return data

            obs_assim   = get_int('assim_flag')
            obs_species = get_int('species')
            obs_tilenum = get_int('tilenum')
            obs_lon     = get_float('lon', mask_nodata=False)
            obs_lat     = get_float('lat', mask_nodata=False)
            obs_obs     = get_float('obs', mask_nodata=False)
            obs_obsvar  = get_float('obsvar')
            obs_fcst    = get_float('fcst')
            obs_fcstvar = get_float('fcstvar')
            obs_ana     = get_float('ana')
            obs_anavar  = get_float('anavar')

    else:
        print(f"file does not exist: {fname}")

    return {'date_time'  : date_time,
            'obs_assim'  : obs_assim,
            'obs_species': obs_species,
            'obs_tilenum': obs_tilenum,
            'obs_lon'    : obs_lon,
            'obs_lat'    : obs_lat,
            'obs_obs'    : obs_obs,
            'obs_obsvar' : obs_obsvar,
            'obs_fcst'   : obs_fcst,
            'obs_fcstvar': obs_fcstvar,
            'obs_ana'    : obs_ana,
            'obs_anavar' : obs_anavar}


# ----------------------------------------------------------------------------
#
# reader for GEOSldas catparam file (binary)

def read_catparam(fname, N_tile, isLDASsa=0):
    """
    Python reader for binary files with Catchment model parameters.
    
    Parameters:
    -----------
    fname : str
        Path to binary Catchment parameter file
    N_tile : int
        Number of tiles; obtain from tilecoord file using read_tilecoord()
    isLDASsa : int, optional
        Flag for LDASsa (1) vs GEOSldas (0) format (default: 0)
        
    Returns:
    --------
    cat_param : dict
        Dictionary containing catchment parameters
    cat_param_units : dict
        Dictionary containing units for each parameter
        
    Notes:
    ------
    Parameter "vegcls" is land cover type:
      vegcls = 1:  broadleaf evergreen trees
      vegcls = 2:  broadleaf deciduous trees
      vegcls = 3:  needleleaf trees
      vegcls = 4:  grassland
      vegcls = 5:  broadleaf shrubs
      vegcls = 6:  dwarf trees
      vegcls = 7:  bare soil
      vegcls = 8:  desert soil
    
    Parameters "soilcls30" and "soilcls100" are 0-30 cm and 0-100 cm soil class.
    
    reichle,  2 Jun 2006
    reichle, 16 Jul 2010 - added vegcls lookup table
    reichle, 28 Oct 2010 - added soilcls*
                         - changed cat_param structure from "vector of
                           structures" to "structure of vectors"
    reichle,  1 Apr 2015 - added new soil parameter fields (file_format==3)
                         - added cat_param_units
    reichle, 28 Jul 2022 - cleaned up LDASsa/GEOSldas switch for commit into GEOSldas repo

    Translated into python from matlab (24 Aug 2026)

    """
    
    # For backward compatibility, back out number of parameters in file
    # from file size:
    # file size = N_param * (N_tile + 2) * bytes_per_datapoint
    
    file_size = os.path.getsize(fname)
    
    if isLDASsa != 0:
        machfmt = '>'  # big-endian, LDASsa
    else:
        machfmt = '<'  # little-endian, GEOSldas
    
    N_param = file_size // ((N_tile + 2) * 4)
    
    if N_param == 40:
        file_format = 1
        if isLDASsa != 0:
            int_columns = [18]  # vegcls (1-indexed in MATLAB, 0-indexed in Python: 17)
        else:
            int_columns = []  # GEOSldas files contain only real*4 numbers
    
    elif N_param in [42, 51, 52]:
        file_format = 2
        if isLDASsa != 0:
            int_columns = [18, 19, 20]  # vegcls, soilcls30, soilcls100
        else:
            int_columns = []  # GEOSldas files contain only real*4 numbers
    
    else:
        raise ValueError('read_catparam: something wrong with file size or format')
    
    print(f'read_catparam: expecting {N_param} parameters in file with file_format {file_format}')
    
    # ----------------------------------------------------------------
    
    int_dtype = np.int32      # precision of fortran tag
    float_dtype = np.float32  # precision of data in input file
    
    print(f'read_catparam: reading from {fname}')
    
    # Adjust int_columns to 0-indexed for Python
    int_columns_py = [i - 1 for i in int_columns]
    
    tmp_data = np.zeros((N_param, N_tile), dtype=np.float32)
    
    with open(fname, 'rb') as ifp:
        
        for i in range(N_param):
            
            # Read fortran tag
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            
            # Read data
            if i in int_columns_py:
                tmp = np.fromfile(ifp, dtype=f'{machfmt}i4', count=N_tile)
            else:
                tmp = np.fromfile(ifp, dtype=f'{machfmt}f4', count=N_tile)
            
            # Read trailing fortran tag
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            
            tmp_data[i, :] = tmp
    
    print('read_catparam: assembling structure array')
    
    cat_param = {}
    cat_param_units = {}
    
    if file_format == 1:
        
        cat_param['dpth'] = tmp_data[0, :]
        cat_param_units['dpth'] = '[mm]'
        
        cat_param['dzsf'] = tmp_data[1, :]
        cat_param_units['dzsf'] = '[mm]'
        cat_param['dzrz'] = tmp_data[2, :]
        cat_param_units['dzrz'] = '[mm]'
        cat_param['dzpr'] = tmp_data[3, :]
        cat_param_units['dzpr'] = '[mm]'
        
        cat_param['dzgt'] = np.column_stack([
            tmp_data[4, :],
            tmp_data[5, :],
            tmp_data[6, :],
            tmp_data[7, :],
            tmp_data[8, :],
            tmp_data[9, :]
        ])
        cat_param_units['dzgt'] = '[m]'
        
        cat_param['poros'] = tmp_data[10, :]
        cat_param_units['poros'] = '[m3 m-3]'
        cat_param['cond'] = tmp_data[11, :]
        cat_param_units['cond'] = '[m s-1]'
        cat_param['psis'] = tmp_data[12, :]
        cat_param_units['psis'] = '[m H2O]'
        cat_param['bee'] = tmp_data[13, :]
        cat_param_units['bee'] = '[-]'
        
        cat_param['wpwet'] = tmp_data[14, :]
        cat_param_units['wpwet'] = '[-]'
        
        cat_param['gnu'] = tmp_data[15, :]
        cat_param_units['gnu'] = '[m-1]'
        
        cat_param['vgwmax'] = tmp_data[16, :]
        cat_param_units['vgwmax'] = '[kg m-2]'
        
        cat_param['vegcls'] = tmp_data[17, :]
        cat_param_units['vegcls'] = '[-]'
        
        cat_param['bf1'] = tmp_data[18, :]
        cat_param_units['bf1'] = '[kg m-4]'
        cat_param['bf2'] = tmp_data[19, :]
        cat_param_units['bf2'] = '[m]'
        cat_param['bf3'] = tmp_data[20, :]
        cat_param_units['bf3'] = '[log(m)]'
        cat_param['cdcr1'] = tmp_data[21, :]
        cat_param_units['cdcr1'] = '[kg m-2]'
        cat_param['cdcr2'] = tmp_data[22, :]
        cat_param_units['cdcr2'] = '[kg m-2]'
        cat_param['ars1'] = tmp_data[23, :]
        cat_param_units['ars1'] = '[m2 kg-1]'
        cat_param['ars2'] = tmp_data[24, :]
        cat_param_units['ars2'] = '[m2 kg-1]'
        cat_param['ars3'] = tmp_data[25, :]
        cat_param_units['ars3'] = '[m4 kg-2]'
        cat_param['ara1'] = tmp_data[26, :]
        cat_param_units['ara1'] = '[m2 kg-1]'
        cat_param['ara2'] = tmp_data[27, :]
        cat_param_units['ara2'] = '[-]'
        cat_param['ara3'] = tmp_data[28, :]
        cat_param_units['ara3'] = '[m2 kg-1]'
        cat_param['ara4'] = tmp_data[29, :]
        cat_param_units['ara4'] = '[-]'
        cat_param['arw1'] = tmp_data[30, :]
        cat_param_units['arw1'] = '[m2 kg-1]'
        cat_param['arw2'] = tmp_data[31, :]
        cat_param_units['arw2'] = '[m2 kg-1]'
        cat_param['arw3'] = tmp_data[32, :]
        cat_param_units['arw3'] = '[m4 kg-2]'
        cat_param['arw4'] = tmp_data[33, :]
        cat_param_units['arw4'] = '[-]'
        cat_param['tsa1'] = tmp_data[34, :]
        cat_param_units['tsa1'] = '[-]'
        cat_param['tsa2'] = tmp_data[35, :]
        cat_param_units['tsa2'] = '[-]'
        cat_param['tsb1'] = tmp_data[36, :]
        cat_param_units['tsb1'] = '[-]'
        cat_param['tsb2'] = tmp_data[37, :]
        cat_param_units['tsb2'] = '[-]'
        cat_param['atau'] = tmp_data[38, :]
        cat_param_units['atau'] = '[-]'
        cat_param['btau'] = tmp_data[39, :]
        cat_param_units['btau'] = '[-]'
    
    elif file_format == 2:
        
        cat_param['dpth'] = tmp_data[0, :]
        cat_param_units['dpth'] = '[mm]'
        
        cat_param['dzsf'] = tmp_data[1, :]
        cat_param_units['dzsf'] = '[mm]'
        cat_param['dzrz'] = tmp_data[2, :]
        cat_param_units['dzrz'] = '[mm]'
        cat_param['dzpr'] = tmp_data[3, :]
        cat_param_units['dzpr'] = '[mm]'
        
        cat_param['dzgt'] = np.column_stack([
            tmp_data[4, :],
            tmp_data[5, :],
            tmp_data[6, :],
            tmp_data[7, :],
            tmp_data[8, :],
            tmp_data[9, :]
        ])
        cat_param_units['dzgt'] = '[m]'
        
        cat_param['poros'] = tmp_data[10, :]
        cat_param_units['poros'] = '[m3 m-3]'
        cat_param['cond'] = tmp_data[11, :]
        cat_param_units['cond'] = '[m s-1]'
        cat_param['psis'] = tmp_data[12, :]
        cat_param_units['psis'] = '[m H2O]'
        cat_param['bee'] = tmp_data[13, :]
        cat_param_units['bee'] = '[-]'
        
        cat_param['wpwet'] = tmp_data[14, :]
        cat_param_units['wpwet'] = '[-]'
        
        cat_param['gnu'] = tmp_data[15, :]
        cat_param_units['gnu'] = '[m-1]'
        
        cat_param['vgwmax'] = tmp_data[16, :]
        cat_param_units['vgwmax'] = '[kg m-2]'
        
        cat_param['vegcls'] = tmp_data[17, :]
        cat_param_units['vegcls'] = '[-]'
        cat_param['soilcls30'] = tmp_data[18, :]
        cat_param_units['soilcls30'] = '[-]'
        cat_param['soilcls100'] = tmp_data[19, :]
        cat_param_units['soilcls100'] = '[-]'
        
        cat_param['bf1'] = tmp_data[20, :]
        cat_param_units['bf1'] = '[kg m-4]'
        cat_param['bf2'] = tmp_data[21, :]
        cat_param_units['bf2'] = '[m]'
        cat_param['bf3'] = tmp_data[22, :]
        cat_param_units['bf3'] = '[log(m)]'
        cat_param['cdcr1'] = tmp_data[23, :]
        cat_param_units['cdcr1'] = '[kg m-2]'
        cat_param['cdcr2'] = tmp_data[24, :]
        cat_param_units['cdcr2'] = '[kg m-2]'
        cat_param['ars1'] = tmp_data[25, :]
        cat_param_units['ars1'] = '[m2 kg-1]'
        cat_param['ars2'] = tmp_data[26, :]
        cat_param_units['ars2'] = '[m2 kg-1]'
        cat_param['ars3'] = tmp_data[27, :]
        cat_param_units['ars3'] = '[m4 kg-2]'
        cat_param['ara1'] = tmp_data[28, :]
        cat_param_units['ara1'] = '[m2 kg-1]'
        cat_param['ara2'] = tmp_data[29, :]
        cat_param_units['ara2'] = '[-]'
        cat_param['ara3'] = tmp_data[30, :]
        cat_param_units['ara3'] = '[m2 kg-1]'
        cat_param['ara4'] = tmp_data[31, :]
        cat_param_units['ara4'] = '[-]'
        cat_param['arw1'] = tmp_data[32, :]
        cat_param_units['arw1'] = '[m2 kg-1]'
        cat_param['arw2'] = tmp_data[33, :]
        cat_param_units['arw2'] = '[m2 kg-1]'
        cat_param['arw3'] = tmp_data[34, :]
        cat_param_units['arw3'] = '[m4 kg-2]'
        cat_param['arw4'] = tmp_data[35, :]
        cat_param_units['arw4'] = '[-]'
        cat_param['tsa1'] = tmp_data[36, :]
        cat_param_units['tsa1'] = '[-]'
        cat_param['tsa2'] = tmp_data[37, :]
        cat_param_units['tsa2'] = '[-]'
        cat_param['tsb1'] = tmp_data[38, :]
        cat_param_units['tsb1'] = '[-]'
        cat_param['tsb2'] = tmp_data[39, :]
        cat_param_units['tsb2'] = '[-]'
        cat_param['atau'] = tmp_data[40, :]
        cat_param_units['atau'] = '[-]'
        cat_param['btau'] = tmp_data[41, :]
        cat_param_units['btau'] = '[-]'
        
        if N_param in [51, 52]:
            
            cat_param['gravel30'] = tmp_data[42, :]
            cat_param_units['gravel30'] = '[%vol]'
            cat_param['orgC30'] = tmp_data[43, :]
            cat_param_units['orgC30'] = '[%weight]'
            cat_param['orgC'] = tmp_data[44, :]
            cat_param_units['orgC'] = '[%weight]'
            cat_param['sand30'] = tmp_data[45, :]
            cat_param_units['sand30'] = '[%weight]'
            cat_param['clay30'] = tmp_data[46, :]
            cat_param_units['clay30'] = '[%weight]'
            cat_param['sand'] = tmp_data[47, :]
            cat_param_units['sand'] = '[%weight]'
            cat_param['clay'] = tmp_data[48, :]
            cat_param_units['clay'] = '[%weight]'
            cat_param['wpwet30'] = tmp_data[49, :]
            cat_param_units['wpwet30'] = '[-]'
            cat_param['poros30'] = tmp_data[50, :]
            cat_param_units['poros30'] = '[m3 m-3]'
        
        if N_param == 52:
            
            cat_param['veghght'] = tmp_data[51, :]
            cat_param_units['veghght'] = '[m]'
    
    else:
        raise ValueError('read_catparam: something wrong with file size or format')
    
    return cat_param, cat_param_units



# ================ EOF =================================================
