#!/usr/bin/env python3

"""
Sample script for plotting maps of long-term data assimilation diagnostics.
Plots Nobs-weighted avg of each metric across all species.
Requires saved files with monthly sums (see Get_ObsFcstAna_sums.py).
Stats of *normalized* OmFs are approximated!

Usage:
    1. Edit "user_config.py" with experiments information.
    2. Run this script as follows (on Discover):
    $ module load python/GEOSpyD
    $ ./Plot_stats_maps.py
"""

import sys;       sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")
import os

import numpy             as np
import matplotlib.pyplot as plt

from datetime               import timedelta

from netCDF4                import Dataset

from tile_to_latlon         import tile_to_latlon, get_grid_resolution
from plot                   import plotMap
from EASEv2                 import EASEv2_ind2latlon

from postproc_ObsFcstAna    import postproc_ObsFcstAna
from user_config            import get_config

def Main_OmF_maps():
    
    # Plot maps of long-term temporal stats
    
    config     = get_config()                 # edit user-defined inputs in user_config.py
    
    exp_list   = config['exp_list']
    start_time = config['start_time']
    end_time   = config['end_time']
    sum_path   = config['sum_path']
    out_path   = config['out_path']

    # min total number of individual O-Fs required for stats at a location (across all months in period)
    
    Nmin = 20
    
    # ------------------------------------------------------------------------------------
    # Read or compute stats.  Each field has dimension [N_tile, N_species] 
    # ------------------------------------------------------------------------------------

    # File name for long-term temporal stats
    stats_file  = out_path + 'temporal_stats_'+exp_list[0]['exptag']+'_'+ start_time.strftime('%Y%m%d')+'_'+ \
        (end_time+timedelta(days=-1)).strftime('%Y%m%d')+'.nc4'

    # Read or compute long-term temporal stats.  Each field has the dimension [N_tile, N_species].
    if os.path.isfile(stats_file):

        print('reading stats nc4 file '+ stats_file)
        stats = {}
        with Dataset(stats_file,'r') as nc:
            for key, value in nc.variables.items():
                stats[key] = value[:].filled(np.nan)

    else:
        # Initialize the postprocessing object and get information
        postproc_obj = postproc_ObsFcstAna(exp_list, start_time, end_time, sum_path=sum_path)
        stats = postproc_obj.calc_temporal_stats_from_sums(write_to_nc=True, fout_stats=stats_file)
    
    N_data        = stats['N_data']
    OmF_mean      = stats['OmF_mean']
    OmF_stdv      = stats['OmF_stdv']
    OmA_mean      = stats['OmA_mean']
    OmA_stdv      = stats['OmA_stdv']
    OmF_norm_mean = stats['OmF_norm_mean']
    OmF_norm_stdv = stats['OmF_norm_stdv']

    # Screen stats with Nmin (except for N_data)
    OmF_mean[     N_data < Nmin] = np.nan
    OmF_stdv[     N_data < Nmin] = np.nan
    OmA_mean[     N_data < Nmin] = np.nan
    OmA_stdv[     N_data < Nmin] = np.nan   
    OmF_norm_mean[N_data < Nmin] = np.nan
    OmF_norm_stdv[N_data < Nmin] = np.nan
    
    # Select/combine fields for plotting. The following provides an example to
    # computes average stats across all species.
    
    # Compute Nobs-weighted avg of each metric across all species.
    # Typically used for SMAP Tb_h/h from asc and desc overpasses,
    # or ASCAT soil moisture from Metop-A/B/C.
    # DOES NOT MAKE SENSE IF, SAY, SPECIES HAVE DIFFERENT UNITS!
    Nobs_data     = np.nansum(              N_data, axis=1)
    OmF_mean      = np.nansum(OmF_mean     *N_data, axis=1)/Nobs_data
    OmF_stdv      = np.nansum(OmF_stdv     *N_data, axis=1)/Nobs_data
    OmF_norm_mean = np.nansum(OmF_norm_mean*N_data, axis=1)/Nobs_data
    OmF_norm_stdv = np.nansum(OmF_norm_stdv*N_data, axis=1)/Nobs_data
    OmA_mean      = np.nansum(OmA_mean     *N_data, axis=1)/Nobs_data
    OmA_stdv      = np.nansum(OmA_stdv     *N_data, axis=1)/Nobs_data

    # Nobs_data need to be at least 1, otherwise unwanted zeros 
    # might get into the computation of global mean etc.
    Nobs_data[Nobs_data == 0] = np.nan  

    # --------------------------------------------------------------------------------
    # Plot stats on regular lat/lon grid, with grid spacing guided by the footprint 
    # scale (field-of-view; FOV) of the  observations.  
    # ---------------------------------------------------------------------------------
    
    # Get geographic info about model grid, model tile space, and obs FOV
    # (assumes that all species have same FOV and FOV_units)

    tc            = exp_list[0]['tilecoord']
    tg            = exp_list[0]['tilegrid_global']
    obs_fov       = exp_list[0]['obsparam'][0]['FOV']
    obs_fov_units = exp_list[0]['obsparam'][0]['FOV_units']

    # Determine spacing of lat/lon grid for plotting (degree) based on obs FOV.
    # Can also set manually if desired.

    my_res = obs_fov*2.      # Note: FOV is a radius, hence the factor of 2.
    
    # If applicable, convert to degree, assuming 100 km is about 1 deg.
    if 'km' in obs_fov_units:
        my_res = my_res/100.

    # If obs_fov=0, reset to match model resolution.
    if obs_fov == 0:
        my_res = np.sqrt( tile_grid['dlat']*tile_grid['dlon'] )

    # pick a resolution from a sample vector of resolutions

    sample_res = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0]

    my_res = min(sample_res, key=lambda x: abs(x - my_res))
    
    print(f'resolution of lat/lon grid for plotting: {my_res} degree')

    # ----------------------------------------------------------------------------------
    # Plot
    # ----------------------------------------------------------------------------------

    exptag = exp_list[0]['exptag']

    fig, axes = plt.subplots(2,2, figsize=(18,10))
    plt.rcParams.update({'font.size':14})
    
    for i in np.arange(2):
        for j in np.arange(2):
            units = '[k]'
            if i == 0 and j == 0:
                tile_data = Nobs_data
                # crange is [cmin, cmax]
                crange =[0, np.ceil((end_time-start_time).days/150)*300]
                colormap = plt.get_cmap('jet',20)
                title_txt = exptag + ' Tb Nobs '    + start_time.strftime('%Y%m')+'_'+end_time.strftime('%Y%m')
                units = '[-]'
            if i == 0 and j ==1:
                tile_data = OmF_mean
                crange =[-3, 3]
                colormap = plt.get_cmap('bwr', 15) 
                title_txt = exptag + ' Tb O-F mean '+ start_time.strftime('%Y%m')+'_'+end_time.strftime('%Y%m')
            if i == 1 and j == 0:
                tile_data = OmF_stdv
                crange =[0, 15]
                colormap = plt.get_cmap ('jet',15)
                title_txt = exptag + ' Tb O-F stdv '+ start_time.strftime('%Y%m')+'_'+end_time.strftime('%Y%m')
            if i == 1 and j == 1:
                tile_data = OmF_norm_stdv
                crange =[0, 15]
                colormap = plt.get_cmap ('jet',15)
                title_txt = exptag + ' Tb normalized O-F stdv (approx!) '+ start_time.strftime('%Y%m%d')+'_'+end_time.strftime('%Y%m%d')

            colormap.set_bad(color='0.9') # light grey, 0-black, 1-white
           
            # map tile_data on 2-d regular lat/lon grid for plotting
            grid_data, lat_2d, lon_2d = tile_to_latlongrid(tile_data, tc, my_res)


            # compute spatial avg of metric (directly from tile_data; no weighting)
            mean    = np.nanmean(       tile_data )
            absmean = np.nanmean(np.abs(tile_data))
            
            if 'normalized' in title_txt:
                absmean = np.nanmean(np.abs(tile_data-1.) )

            if 'normalized' in title_txt and 'stdv' in title_txt:
                title_txt = title_txt + '\n' + "avg=%.3f, avg(abs(nstdv-1))=%.3f" % (mean, absmean)+' '+units
            elif 'mean' in title_txt:
                title_txt = title_txt + '\n' + "avg=%.3f, avg(abs)=%.3f"          % (mean, absmean)+' '+units
            else:
                title_txt = title_txt + '\n' + "avg=%.2f"                         % (mean         )+' '+units                
   
            if 'normalized' in title_txt:
                grid_data = np.log10(grid_data)
                crange = [-0.6, 0.45]
            
            mm, cs = plotMap(grid_data, ax =axes[i,j], lat=lat_2d, lon=lon_2d, cRange=crange, \
                        title=title_txt, cmap=colormap, bounding=[-60, 80, -180,180])

    plt.tight_layout()
    # Save figure to file
    fig.savefig(out_path+'Map_OmF_'+ exptag +'_'+start_time.strftime('%Y%m')+'_'+\
                        (end_time+timedelta(days=-1)).strftime('%Y%m')+'.png')
    #plt.show()
    plt.close(fig)

# -----------------------------------------------------------------------------------------------
    
if __name__ == '__main__':
        
    Main_OmF_maps()

# ====================== EOF =========================================================
