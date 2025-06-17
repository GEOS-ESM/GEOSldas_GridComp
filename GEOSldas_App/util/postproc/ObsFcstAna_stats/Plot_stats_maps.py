#!/usr/bin/env python3

"""
Sample script for plotting maps of long-term data assimilation diagnostics.
Plots Nobs-weighted avg of each metric across all species.
Requires saved files with monthly sums (see Get_ObsFcstAna_stat.py).
Stats of *normalized* OmFs are approximated!
"""

import sys;       sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")

import numpy             as np
import matplotlib.pyplot as plt

from datetime               import timedelta

from tile2grid              import tile2grid
from plot                   import plotMap
from EASEv2                 import EASEv2_ind2latlon

def plot_OmF_maps(postproc_obj, stats, fig_path='./'):
    
    start_time    = postproc_obj.start_time
    end_time      = postproc_obj.end_time
    domain        = postproc_obj.domain
    tc            = postproc_obj.tilecoord
    tg            = postproc_obj.tilegrid_global
    
    N_data        = stats['N_data']
    OmF_mean      = stats['OmF_mean']
    OmF_stdv      = stats['OmF_stdv']
    OmA_mean      = stats['OmA_mean']
    OmA_stdv      = stats['OmA_stdv']
    OmF_norm_mean = stats['OmF_norm_mean']
    OmF_norm_stdv = stats['OmF_norm_stdv']

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

    # Plotting
    exptag = postproc_obj.exptag

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

            # Regrid 1d tile_data to 2d grid_data for map plots
            if '_M09_' in domain: # special case  
                grid_data_M09 = np.zeros((1624, 3856)) + np.nan  
                grid_data_M09[tc['j_indg'],tc['i_indg']] = tile_data
                
                # Reshape the data into 4x4 blocks
                reshaped = grid_data_M09.reshape(1624//4, 4, 3856//4, 4)

                # Combine each 4x4 M09 block into a M36 grid
                #if i==0 and j==0:
                #   grid_data = np.sum(reshaped,axis=(1, 3)) 
                #else:
                #   grid_data = np.nanmean(reshaped,axis=(1, 3))

                grid_data = grid_data_M09[1::4, 2::4]

                # NOT area weighted 
                wmean = np.nanmean(grid_data)
                wabsmean = np.nanmean(np.abs(grid_data))
                if 'normalized' in title_txt:
                    wabsmean = np.nanmean(np.abs(grid_data-1.))
                    
                lat_M36, lon_M36 = EASEv2_ind2latlon(np.arange(406), np.arange(964),'M36')
                lon_2d,lat_2d = np.meshgrid(lon_M36,lat_M36)
            else:
                grid_data, lat_2d, lon_2d = tile2grid(tile_data, tc, tg)
                
                # Area weighted mean and mean(abs)
                wmean    =     np.nansum(       tile_data     * tc['area'])/np.nansum(~np.isnan(tile_data)*tc['area'])
                wabsmean =     np.nansum(np.abs(tile_data)    * tc['area'])/np.nansum(~np.isnan(tile_data)*tc['area'])
                if 'normalized' in title_txt:
                    wabsmean = np.nansum(np.abs(tile_data-1.) * tc['area'])/np.nansum(~np.isnan(tile_data)*tc['area'])
                    
            if 'normalized' in title_txt:
                title_txt = title_txt + '\n' + "avg=%.3f, avg(abs(nstdv-1))=%.3f" % (wmean, wabsmean)+' '+units
            elif 'mean' in title_txt:
                title_txt = title_txt + '\n' + "avg=%.3f, avg(abs)=%.3f" % (wmean, wabsmean)+' '+units
            else:
                title_txt = title_txt + '\n' + "avg=%.2f" % (wmean) +' '+units                
         
            if 'normalized' in title_txt:
                grid_data = np.log10(grid_data)
                crange = [-0.6, 0.45]
            
            mm, cs = plotMap(grid_data, ax =axes[i,j], lat=lat_2d, lon=lon_2d, cRange=crange, \
                        title=title_txt, cmap=colormap, bounding=[-60, 80, -180,180])            
            
    plt.tight_layout()
    # Save figure to file
    fig.savefig(fig_path+'Map_OmF_'+ exptag +'_'+start_time.strftime('%Y%m')+'_'+\
                        (end_time+timedelta(days=-1)).strftime('%Y%m')+'.png')
    #plt.show()
    plt.close(fig)

# -----------------------------------------------------------------------------------------------
    
if __name__ == '__main__':

    # Plot maps of long-term temporal stats
    
    from postproc_ObsFcstAna    import postproc_ObsFcstAna
    from user_config            import get_config
    
    config     = get_config()                 # edit user-defined inputs in user_config.py
    
    exp_list   = config['exp_list']
    start_time = config['start_time']
    end_time   = config['end_time']
    sum_path   = config['sum_path']
    out_path   = config['out_path']
     
    # ------------------------------------------------------------------------------------

    # Initialize the postprocessing object
    postproc = postproc_ObsFcstAna(exp_list, start_time, end_time, sum_path=sum_path)

    # File name for long-term temporal stats
    stats_file  = out_path + 'temporal_stats_'+exp_list[0]['exptag']+'_'+ start_time.strftime('%Y%m%d')+'_'+ \
        (end_time+timedelta(days=-1)).strftime('%Y%m%d')+'.nc4'


    # Read or compute long-term temporal stats.  Each field has the dimension [N_tile, N_species].
    if os.path.isfile(stats_file):

        print('reading stats nc4 file '+ stats_file)
        temporal_stats = {}
        with Dataset(stats_file,'r') as nc:
            for key, value in nc.variables.items():
                temporal_stats[key] = value[:].filled(np.nan)

    else:
        temporal_stats = postproc.calc_temporal_stats_from_sums(write_to_nc=True, fout_stats=stats_file)

    # Plot stats
        
    plot_OmF_maps(postproc, temporal_stats, fig_path=out_path )

# ====================== EOF =========================================================
