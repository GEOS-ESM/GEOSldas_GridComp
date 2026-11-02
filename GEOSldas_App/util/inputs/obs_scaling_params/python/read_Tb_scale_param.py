import numpy as np
import struct

def read_Tb_scale_param(fname):
    """
    Python reader for binary files with GEOSldas Tb obs scaling params used in SMAP L4_SM algorithm
    (see GEOSldas Fortran subroutine scale_obs_Tb_zscore()).
    
    Example file name:  L4SM_v015_L1C_zscore_stats_D_p04.bin
    
    Stores version 015 (v015) scaling parameters (used in L4_SM Version 8) for SMAP Tb
    observations from descending (D) half-orbits for pentad 4 (p04; 16-20 January). 
    
    The L4_SM algorithm runs on the 9-km EASEv2_M09 grid and assimilates SMAP L1C_Tb observations,
    which are on the EASEv2_M36 grid.
    The Tb observations are "scaled" to the climatology of the simulated (36-km) Tbs by subtracting
    the climatological mean of the SMAP Tb observations and adding the climatological mean of the
    corresponding simulated Tbs. The climatological mean values are seasonally varying and stored
    in pentad scaling files.
    
    Parameters:
    -----------
    fname : str
        Path to the binary file
        
    Returns:
    --------
    asc_flag : int
        Orbit direction: 1=ascending (6pm), 0=descending (6am)
    N_sclprm : int
        Number of elements (L1C_TB EASEv2_M36 grid cells) in scaling parameter vectors
    N_angle : int
        Number of Tb incidence angles (=1)
    inc_angle : numpy.ndarray
        Tb incidence angle (=40)
    lon : numpy.ndarray
        Center longitude of 9-km grid cell associated with Tb scaling parameters
    lat : numpy.ndarray
        Center latitude  of 9-km grid cell associated with Tb scaling parameters
    tbh_mean_obs : numpy.ndarray
        Mean of observed  H-pol Tb for pentad and grid cell
    tbh_std_obs : numpy.ndarray
        Stdv of observed  H-pol Tb for pentad and grid cell [not used in L4_SM]
    tbh_mean_mod : numpy.ndarray
        Mean of simulated H-pol Tb for pentad and grid cell
    tbh_std_mod : numpy.ndarray
        Stdv of simulated H-pol Tb for pentad and grid cell [not used in L4_SM]
    tbh_N_data : numpy.ndarray
        Number of data used in computation of H-pol mean and stdv values
    tbv_mean_obs : numpy.ndarray
        Mean of observed  V-pol Tb for pentad and grid cell
    tbv_std_obs : numpy.ndarray
        Stdv of observed  V-pol Tb for pentad and grid cell [not used in L4_SM]
    tbv_mean_mod : numpy.ndarray
        Mean of simulated V-pol Tb for pentad and grid cell
    tbv_std_mod : numpy.ndarray
        Stdv of simulated V-pol Tb for pentad and grid cell [not used in L4_SM]
    tbv_N_data : numpy.ndarray
        Number of data used in computation of V-pol mean and stdv values
    
    - reichle, 24 Aug 2026
    """
    
    int_dtype    = np.int32      # precision of fortran tag and integer data
    float_dtype  = np.float32    # precision of data in input file
    
    nodata_stats = -9999         # no-data-value for stats data
    
    print()
    print(f'reading from {fname}')
    
    # Open file in binary read mode with big-endian byte order
    
    with open(fname, 'rb') as ifp:
        
        # Read header lines:
        
        # Header line 1: asc_flag (orbit direction), N_data_min
        fortran_tag = struct.unpack('>i',  ifp.read(4))[0]
        tmp_data    = struct.unpack('>2i', ifp.read(8))
        fortran_tag = struct.unpack('>i',  ifp.read(4))[0]
        
        asc_flag    = tmp_data[0]
        N_data_min  = tmp_data[1]  # ignore, see "tb[X]_N_data"
        
        print()
        print(f'asc_flag                     : {asc_flag}')
        print()
        print(f'[ignore] N_data_min          : {N_data_min}')
        
        # Header line 2: start year/month/day/hour/min
        fortran_tag = struct.unpack('>i',  ifp.read(4))[0]
        tmp_data    = struct.unpack('>5i', ifp.read(20))
        fortran_tag = struct.unpack('>i',  ifp.read(4))[0]
        
        start_time = {
            'year':  tmp_data[0],
            'month': tmp_data[1],
            'day':   tmp_data[2],
            'hour':  tmp_data[3],
            'min':   tmp_data[4]
        }
        
        print(f'[ignore] start_time.year     : {start_time["year"]}')
        print(f'[ignore] start_time.month    : {start_time["month"]}')
        print(f'[ignore] start_time.day      : {start_time["day"]}')
        print(f'[ignore] start_time.hour     : {start_time["hour"]}')
        print(f'[ignore] start_time.min      : {start_time["min"]}')
        
        # Header line 3: end year/month/day/hour/min
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tmp_data = struct.unpack('>5i', ifp.read(20))
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        end_time = {
            'year': tmp_data[0],
            'month': tmp_data[1],
            'day': tmp_data[2],
            'hour': tmp_data[3],
            'min': tmp_data[4]
        }
        
        print(f'[ignore] end_time.year       : {end_time["year"]}')
        print(f'[ignore] end_time.month      : {end_time["month"]}')
        print(f'[ignore] end_time.day        : {end_time["day"]}')
        print(f'[ignore] end_time.hour       : {end_time["hour"]}')
        print(f'[ignore] end_time.min        : {end_time["min"]}')
        print()
        
        # Header line 4: N_sclprm, N_angle, N_tile
        fortran_tag = struct.unpack('>i',  ifp.read(4))[0]
        tmp_data    = struct.unpack('>3i', ifp.read(12))
        fortran_tag = struct.unpack('>i',  ifp.read(4))[0]
        
        N_sclprm = tmp_data[0]  # = N_tile for tile-space scaling params
        N_angle  = tmp_data[1]
        N_tile   = tmp_data[2]
        
        print(f'N_sclprm                     : {N_sclprm}')
        print(f'N_angle                      : {N_angle}')
        print()
        print(f'[ignore] N_tile              : {N_tile}')
        print()
        
        if N_angle != 1:
            raise ValueError('reader only works for N_angle=1')
        
        # Header line 5: incidence angle(s)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        inc_angle = np.fromfile(ifp, dtype='>f4', count=N_angle)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        print(f'inc_angle(s)                 : {inc_angle}')
        
        # Coord info:
        
        # Records 6-7: longitude and latitude
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        lon = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        lat = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        print(f'min max lon                  : {[np.min(lon), np.max(lon)]}')
        print(f'min max lat                  : {[np.min(lat), np.max(lat)]}')
        
        # Record 8: tile ID (EASEv2_M09 tile ID associated with Tb scaling parameter)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tile_id = np.fromfile(ifp, dtype='>i4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        print()
        print(f'[ignore] min max tile_id     : {[np.min(tile_id), np.max(tile_id)]}')
        print()
        
        # Records 9-13: TbH mean_obs, std_obs, mean_mod, std_mod, N_data
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tbh_mean_obs = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tbh_std_obs = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tbh_mean_mod = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tbh_std_mod = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tbh_N_data = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        # Records 14-18: TbV mean_obs, std_obs, mean_mod, std_mod, N_data
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tbv_mean_obs = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tbv_std_obs = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tbv_mean_mod = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tbv_std_mod = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tbv_N_data = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        # Records 19-22: unused data of length N_sclprm
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tmprecord19 = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tmprecord20 = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tmprecord21 = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
        tmprecord22 = np.fromfile(ifp, dtype='>f4', count=N_sclprm)
        fortran_tag = struct.unpack('>i', ifp.read(4))[0]
    
    # Replace no-data-values with NaNs
    tbh_mean_obs[tbh_mean_obs == nodata_stats] = np.nan
    tbh_std_obs[tbh_std_obs == nodata_stats] = np.nan
    tbh_mean_mod[tbh_mean_mod == nodata_stats] = np.nan
    tbh_std_mod[tbh_std_mod == nodata_stats] = np.nan
    tbh_N_data[tbh_N_data == nodata_stats] = np.nan
    
    tbv_mean_obs[tbv_mean_obs == nodata_stats] = np.nan
    tbv_std_obs[tbv_std_obs == nodata_stats] = np.nan
    tbv_mean_mod[tbv_mean_mod == nodata_stats] = np.nan
    tbv_std_mod[tbv_std_mod == nodata_stats] = np.nan
    tbv_N_data[tbv_N_data == nodata_stats] = np.nan
    
    tmprecord19[tmprecord19 == nodata_stats] = np.nan
    tmprecord20[tmprecord20 == nodata_stats] = np.nan
    tmprecord21[tmprecord21 == nodata_stats] = np.nan
    tmprecord22[tmprecord22 == nodata_stats] = np.nan
    
    print(f'min max tbh_mean_obs         : {[np.nanmin(tbh_mean_obs), np.nanmax(tbh_mean_obs)]}')
    print(f'min max tbh_std_obs          : {[np.nanmin(tbh_std_obs),  np.nanmax(tbh_std_obs )]}')
    print(f'min max tbh_mean_mod         : {[np.nanmin(tbh_mean_mod), np.nanmax(tbh_mean_mod)]}')
    print(f'min max tbh_std_mod          : {[np.nanmin(tbh_std_mod),  np.nanmax(tbh_std_mod )]}')
    print(f'min max tbh_N_data           : {[np.nanmin(tbh_N_data),   np.nanmax(tbh_N_data  )]}')
    print()
    
    print(f'min max tbv_mean_obs         : {[np.nanmin(tbv_mean_obs), np.nanmax(tbv_mean_obs)]}')
    print(f'min max tbv_std_obs          : {[np.nanmin(tbv_std_obs),  np.nanmax(tbv_std_obs )]}')
    print(f'min max tbv_mean_mod         : {[np.nanmin(tbv_mean_mod), np.nanmax(tbv_mean_mod)]}')
    print(f'min max tbv_std_mod          : {[np.nanmin(tbv_std_mod),  np.nanmax(tbv_std_mod )]}')
    print(f'min max tbv_N_data           : {[np.nanmin(tbv_N_data),   np.nanmax(tbv_N_data  )]}')
    print()
    
    print(f'[ignore] min max tmprecord19 : {[np.nanmin(tmprecord19),  np.nanmax(tmprecord19) ]}')
    print(f'[ignore] min max tmprecord20 : {[np.nanmin(tmprecord20),  np.nanmax(tmprecord20) ]}')
    print(f'[ignore] min max tmprecord21 : {[np.nanmin(tmprecord21),  np.nanmax(tmprecord21) ]}')
    print(f'[ignore] min max tmprecord22 : {[np.nanmin(tmprecord22),  np.nanmax(tmprecord22) ]}')
    print()
    
    return (asc_flag, N_sclprm, N_angle, inc_angle,
            lon, lat,
            tbh_mean_obs, tbh_std_obs, tbh_mean_mod, tbh_std_mod, tbh_N_data,
            tbv_mean_obs, tbv_std_obs, tbv_mean_mod, tbv_std_mod, tbv_N_data)

# --------------- EOF -------------------------------------------------------------------
