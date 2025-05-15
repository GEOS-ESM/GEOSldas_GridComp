#!/usr/bin/env python3

"""
Sample script for computing/plotting data assimilation (DA) diagnostics from 
pre-saved monthly sums.

"""

import sys;       sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")
import os

import numpy             as np
import matplotlib.pyplot as plt

from datetime               import datetime, timedelta
from dateutil.relativedelta import relativedelta
from netCDF4                import Dataset, num2date
from mpl_toolkits.basemap   import Basemap

from read_GEOSldas          import read_tilecoord, read_obs_param
from tile2grid             import tile2grid
from plot                   import plotMap
from EASEv2                 import EASEv2_ind2latlon

from postproc_ObsFcstAna    import postproc_ObsFcstAna

# Uncomment if to run the script in the background to see the standard output while running 
# import io
#sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
#sys.stderr = io.TextIOWrapper(open(sys.stderr.fileno(), 'wb', 0), write_through=True)

# User provided time range for processing (Year, Month, Day)
start_time = datetime(2015,4,1)
end_time   = datetime(2016,4,1)

# User provided output directory 
out_path = '/discover/nobackup/qliu/SMAP_diag/'
os.makedirs(out_path, exist_ok=True)

# -------------------------------- Experiments Information -----------------------------------------
# Supports single experiment or multiple experiments.
# All experiments must have identical tilecoords and number/order of observation species.
# If the default "species" number/order do not match, need to set the *optional*
#   "select_species" key to get a match, i.e. same species sequences.
# This capability is required to enable calculating OmF/OmA statistics for one experiment
#   using observations from another experiment. See note below.

exp_main = { 'expdir'       : '/discover/nobackup/projects/gmao/merra/iau/merra_land/SMAP_runs/SMAP_Nature_v11/',
          'expid'        : 'DAv8_SMOSSMAP',
          'exptag'       : 'DAMulti_SMAP', 
          'domain'       : 'SMAP_EASEv2_M36_GLOBAL',
          'species_list' : [5,6,7,8] }

exp_sup1 = { 'expdir'       : '/discover/nobackup/projects/gmao/merra/iau/merra_land/SMAP_runs/SMAP_Nature_v11/',
          'expid'        : 'DAv8_M36',
          'exptag'       : 'DASMAP_SMAP', 
          'domain'       : 'SMAP_EASEv2_M36_GLOBAL',
          'species_list' : [1,2,3,4]  }

# User provided monthly sum files directory
monthly_sums_path = '/discover/nobackup/qliu/SMAP_diag/' +exp_main['expid']+ \
              '/output/'+exp_main['domain']+'/ana/ens_avg/'

# Uses forecasts/analyses from first experiment in list.
# Observations from experiment specified by 'obs_from' index.
# The mostly likely use case for this is that _scaled_ observations from a DA experiment
#   are used to compute OmF etc diagnostics for a corresponding open loop experiment.

exp_list = [exp_main, exp_sup1]
obs_from = 1                            # obs is from "exp_sup1" (0-based indexing)
if obs_from >= len(exp_list):
    raise ValueError('Invalid "obs_from" value')

# Add tilecoord and obs_param information to each experiment
for exp in exp_list:
    expdir   = exp['expdir']
    expid    = exp['expid']
    domain   = exp['domain']
    fop      = expdir+expid+'/output/'+domain+'/rc_out/Y2015/M04/'+expid+'.ldas_obsparam.20150401_0000z.txt'
    obsparam = read_obs_param(fop)

    # get the species list and default to list of all species if doesn't exist 
    species_list = exp.get('species_list',[int(obsparam[i]['species']) for i in range(len(obsparam))])
    
    # reorder obsparam to match across experiments
    obsparam_new = []
    for i in range(len(obsparam)):
        if int(obsparam[i]['species']) in species_list:
               obsparam_new.append(obsparam[i])              
    obsparam = obsparam_new
    
    ftc = expdir+expid+'/output/'+ domain+'/rc_out/'+ expid+'.ldas_tilecoord.bin'
    tc = read_tilecoord(ftc)

    exp.update({'tilecoord':tc,'obsparam':obsparam})

if len(exp_list) >1 :
    stats_file  = out_path + 'tmp_stats_'+exp_list[0]['exptag']+'_obsfrom_'+ \
                  exp_list[obs_from]['exptag']+'_'+start_time.strftime('%Y%m%d')+'_'+ \
                  end_time.strftime('%Y%m%d')+'.nc4'
else:
    stats_file  = out_path + 'tmp_stats_'+exp_list[0]['exptag']+'_'+ start_time.strftime('%Y%m%d')+'_'+ \
                  end_time.strftime('%Y%m%d')+'.nc4'

#  =========================================================================
#  Postprocess raw ObsFcstAna output data into monthly sums for simpler and faster postprocessing;
#  computes mean, variance from monthly sums that can be used to compute DA diagnostics directly

if not os.path.isfile(stats_file):
    # Initialize the postprocessing object
    postproc = postproc_ObsFcstAna(exp_list, start_time, end_time, obs_from=obs_from)
    # Step 1: Compute and save monthly sums 
    postproc.save_monthly_sum(monthly_sums_path)
    # Step 2: Compute statistics from monthly sums, option to save result to file
    stats = postproc.calculate_stats_from_sums(mo_path=monthly_sums_path, write_to_nc=True, filename=stats_file)
else:
    print('reading stats nc4 file '+stats_file)
    stats = {}
    with Dataset(stats_file,'r') as nc:
        for key, value in nc.variables.items():
            stats[key] = value[:].filled(np.nan)
# 
#  ==========================================================================

# Sample of final compuation of selected diagnostic metrics 
 
Nmin = 20

# Then computer metrics of O-F, O-A, etc. based on above computed
N_data = stats['N_data']
O_mean = stats['obs_mean']
# mean(x-y) = E[x] - E[y]   
OmF_mean = stats['obs_mean'] - stats['fcst_mean']
OmA_mean = stats['obs_mean'] - stats['ana_mean']
# var(x-y) = var(x) + var(y) - 2cov(x,y)
# cov(x,y) = E[xy] - E[x]E[y]
OmF_stdv  = np.sqrt(stats['obs_variance'] + stats['fcst_variance'] - \
                    2 * (stats['oxf_mean'] - stats['obs_mean']*stats['fcst_mean']))
                    
OmA_stdv  = np.sqrt(stats['obs_variance'] + stats['ana_variance'] - \
                    2 * (stats['oxa_mean'] - stats['obs_mean']*stats['ana_mean']))

 # "fcstvar" is assumed constant here for convenience. Modify if necessary
OmF_norm_mean = OmF_mean / np.sqrt(stats['obsvar_mean'] + stats['fcstvar_mean']) 
OmF_norm_stdv = np.sqrt(OmF_stdv**2 / (stats['obsvar_mean'] + stats['fcstvar_mean']) )
  
# Mask out data points with insufficent observations using the Nmin threshold
# Do NOT apply to N_data
OmF_mean[     N_data < Nmin] = np.nan
OmF_stdv[     N_data < Nmin] = np.nan
OmF_norm_mean[N_data < Nmin] = np.nan
OmF_norm_stdv[N_data < Nmin] = np.nan
OmA_mean[     N_data < Nmin] = np.nan
OmA_stdv[     N_data < Nmin] = np.nan
N_data[       N_data < Nmin] = 0

# Combine metrics of individual species using weighted averaging
Nobs_data     = np.nansum(              N_data, axis=1)
OmF_mean      = np.nansum(OmF_mean     *N_data, axis=1)/Nobs_data
OmF_stdv      = np.nansum(OmF_stdv     *N_data, axis=1)/Nobs_data
OmF_norm_mean = np.nansum(OmF_norm_mean*N_data, axis=1)/Nobs_data
OmF_norm_stdv = np.nansum(OmF_norm_stdv*N_data, axis=1)/Nobs_data
OmA_mean      = np.nansum(OmA_mean     *N_data, axis=1)/Nobs_data
OmA_stdv      = np.nansum(OmA_stdv     *N_data, axis=1)/Nobs_data

# Plotting
expid = exp_list[0]['exptag']

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
            title_txt = expid + ' Tb Nobs '+ start_time.strftime('%Y%m')+'_'+end_time.strftime('%Y%m')
            units = '[-]'
        if i == 0 and j ==1:
            tile_data = OmF_mean
            crange =[-3, 3]
            colormap = plt.get_cmap('bwr', 15) 
            title_txt = expid + ' Tb O-F mean '+ start_time.strftime('%Y%m')+'_'+end_time.strftime('%Y%m')
        if i == 1 and j == 0:
            tile_data = OmF_stdv
            crange =[0, 15]
            colormap = plt.get_cmap ('jet',15)
            title_txt = expid + ' Tb O-F stdv '+ start_time.strftime('%Y%m')+'_'+end_time.strftime('%Y%m')
        if i == 1 and j == 1:
            tile_data = OmF_norm_stdv
            crange =[0, 15]
            colormap = plt.get_cmap ('jet',15)
            title_txt = expid + ' Tb normalized O-F stdv '+ start_time.strftime('%Y%m%d')+'_'+end_time.strftime('%Y%m%d')

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
            grid_data, uy, ux = tile2grid(tile_data, lat_1d = tc['com_lat'], lon_1d = tc['com_lon'])
            lon_2d,lat_2d = np.meshgrid(ux, uy)

            # Aear weighted mean and mean(abs)
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
fig.savefig(out_path+'Map_OmF_'+ expid +'_'+start_time.strftime('%Y%m')+'_'+\
                    end_time.strftime('%Y%m')+'.png')
plt.show()
plt.close(fig)


# ====================== EOF =========================================================
