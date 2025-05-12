
# collection of readers for GEOSldas output files

import numpy as np
import struct
import os
import numpy as np

# ----------------------------------------------------------------------------

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

def read_tilecoord(fname):
    int_precision = 'i'
    float_precision = 'f'

    # Determine endianness
    machfmt = '<'  # '>' for big-endian, '<' for little-endian

    print(f"reading from {fname}")

    tile_coord = {}

    with open(fname, 'rb') as ifp:
        fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
        tile_coord['N_tile'] = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
        fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

        Nt = tile_coord['N_tile']

        fields = ['tile_id', 'typ', 'pfaf', 'com_lon', 'com_lat', 'min_lon', 'max_lon',
                  'min_lat', 'max_lat', 'i_indg', 'j_indg', 'frac_cell', 'frac_pfaf',
                  'area', 'elev']

        for field in fields:
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            dtype = int_precision if field in ['tile_id', 'typ', 'pfaf', 'i_indg', 'j_indg'] else float_precision
            tile_coord[field] = np.frombuffer(ifp.read(Nt * 4), dtype=f'{machfmt}{dtype}')
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

    print("done reading file")
    return tile_coord

# ----------------------------------------------------------------------------

def read_ObsFcstAna(fname, isLDASsa=False):

    # Initialize outputs
    nodata = -9999

    date_time = {
        'year': nodata,
        'month': nodata,
        'day': nodata,
        'hour': nodata,
        'min': nodata,
        'sec': nodata,
        'dofyr': nodata,
        'pentad': nodata
    }

    obs_assim = []
    obs_species = []
    obs_tilenum = []
    obs_lon = []
    obs_lat = []
    obs_obs = []
    obs_obsvar = []
    obs_fcst = []
    obs_fcstvar = []
    obs_ana = []
    obs_anavar = []

    # Determine endianness
    machfmt = '>' if isLDASsa else '<'  # '>' for big-endian, '<' for little-endian

    if os.path.exists(fname):
        print(f"reading from {fname}")

        with open(fname, 'rb') as ifp:
            # Read N_obs and time stamp entry
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            N_obs = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            year, month, day, hour, minute, second, dofyr, pentad = struct.unpack(f'{machfmt}8i', ifp.read(32))
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            date_time = {
                'year': year,
                'month': month,
                'day': day,
                'hour': hour,
                'min': minute,
                'sec': second,
                'dofyr': dofyr,
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
            obs_lon = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read latitude
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_lat = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation value
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_obs = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation variance
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_obsvar = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space model forecast value
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_fcst = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space model forecast variance
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_fcstvar = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space analysis value
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_ana = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space analysis variance
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_anavar = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # No-data check
            obs_obsvar[obs_obsvar == nodata] = np.nan
            obs_fcst[obs_fcst == nodata] = np.nan
            obs_fcstvar[obs_fcstvar == nodata] = np.nan
            obs_ana[obs_ana == nodata] = np.nan
            obs_anavar[obs_anavar == nodata] = np.nan

    else:
        print(f"file does not exist: {fname}")

    return {'date_time': date_time, 
            'obs_assim': obs_assim, 
            'obs_species': obs_species, 
            'obs_tilenum': obs_tilenum, 
            'obs_lon': obs_lon, 
            'obs_lat': obs_lat,
            'obs_obs': obs_obs, 
            'obs_obsvar': obs_obsvar, 
            'obs_fcst': obs_fcst, 
            'obs_fcstvar': obs_fcstvar, 
            'obs_ana': obs_ana, 
            'obs_anavar': obs_anavar}

# ================ EOF =================================================
