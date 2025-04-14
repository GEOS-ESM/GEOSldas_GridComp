#!/usr/bin/env python3

"""
This script provids a sample to computes and plots GEOSldas Data Assimilation (DA)
diagnostics maps based on the ObsFcstAna output.

User provides the experiment inforamtion and calls ObsFcstAna_prep to
1. calculate and save the monthly sums of Obs/Fcst/Ana and sum of squared/products
2. calculate mean/stdv of Obs/Fcst/Ana

Afterwards, user can compute the desired DA diagnostic statistics based on the priliminary
results  and create plots. 

To run this script on  Discover:
    $ module load python/GEOSpyD
    $ ./Main_example.py
     or to run in the background,
    $ nohup ./Main_example.py > out.log &

Requirements:
    - Python 3.x
    - Modules: numpy, matplotlib, netCDF4 (included in GEOSpyD)

Author: Q. Liu
Last Modified: Apr., 2025
"""
import sys
sys.path.append('../../shared/python/')
import warnings; warnings.filterwarnings("ignore")
import os

from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from read_GEOSldas import read_tilecoord, read_obs_param
from util import make_folder, array2grid
from plot import plotMap
from easev2 import easev2_ind2latlon

from ObsFcstAna_prep import obsfcstana_prep

# Uncomment if running in the background
# import io
#sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
#sys.stderr = io.TextIOWrapper(open(sys.stderr.fileno(), 'wb', 0), write_through=True)

# -------------------------------- Information of experiment  -----------------------------------------
#  Experiment dictionary, need "expdir", "expid", "domain", "tilecoord", "obsparam"
Vv8010 = { 'expdir' : '/gpfsm/dnb05/projects/p51/SMAP_Nature/',
                    'expid' : 'SPL4SM_Vv8010',
                    'domain':  'SMAP_EASEv2_M09_GLOBAL'}

start_time = datetime(2015,4,1)
end_time = datetime(2021,4,1)

 # Read tilecoord and obsparam for tile and obs species information
expdir = Vv8010['expdir']
expid = Vv8010['expid']
domain = Vv8010['domain']
ftc = expdir+expid+'/output/'+ domain+'/rc_out/'+ expid+'.ldas_tilecoord.bin'
tc = read_tilecoord(ftc)

fop = expdir+expid+'/output/'+domain+'/rc_out/Y2015/M04/'+expid+'.ldas_obsparam.20150401_0000z.txt'
obsparam = read_obs_param(fop)

Vv8010.update({'tilecoord':tc,'obsparam':obsparam})

# Base directory for storing monthly files
# This can be the same as the experiment directory (expdir) or a different location
out_path_mo = '/discover/nobackup/qliu/SMAP_diag/' +Vv8010['expid']+ \
              '/output/'+Vv8010['domain']+'/ana/ens_avg/'

# Output directory for final diagnostic files and plots
out_path = '/discover/nobackup/qliu/SMAP_diag/'
make_folder(out_path)
stats_file  = out_path + 'tmp_stats_Vv8010_'+start_time.strftime('%Y%m%d')+'_'+ \
              end_time.strftime('%Y%m%d')+'.nc4'

#  --------- Preprocess into month data and compute basic statistics --------------------
if not os.path.isfile(stats_file):
    # Initialize
    prep_v8 = obsfcstana_prep(Vv8010, start_time, end_time)
    # Step 1: save monthly sums if files don't exist
    prep_v8.save_monthly_sum(out_path_mo)
    # Step 2: Compute basic statistics based on monthly sums, option to save stats in a nc4 file
    stats = prep_v8.calculate_stats_fromsums(mo_path=out_path_mo, write_to_nc=True, filename=stats_file)
else:
    # Read from previous saved stats file
    print('reading stats nc4 file '+stats_file)
    stats = {}
    with Dataset(stats_file,'r') as nc:
        for key, value in nc.variables.items():
            stats[key] = value[:].filled(np.nan)

#  -------------   Calculate user prefered statistical metrics and make plots ----------------
# Define a minimum threshold for the temporal data points to ensure statistical reliability
# of the computed metrics. 
Nmin = 20

# Then computer metrics of O-F, O-A, etc. based on above computed
N_data = stats['N_data']
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
OmF_mean[N_data < Nmin] = np.nan
OmF_stdv[N_data < Nmin] = np.nan
OmF_norm_mean[N_data < Nmin] = np.nan
OmF_norm_stdv[N_data < Nmin] = np.nan
OmA_mean[N_data < Nmin] = np.nan
OmA_stdv[N_data < Nmin] = np.nan

# Combine metrics of individual species using weighted averaging
OmF_mean = np.nansum(OmF_mean*N_data, axis=1)/np.nansum(N_data,axis=1)
OmF_stdv = np.nansum(OmF_stdv*N_data,axis=1)/np.nansum(N_data,axis=1)
OmF_norm_mean = np.nansum(OmF_norm_mean*N_data, axis=1)/np.nansum(N_data,axis=1)
OmF_norm_stdv = np.nansum(OmF_norm_stdv*N_data,axis=1)/np.nansum(N_data,axis=1)
OmA_mean = np.nansum(OmA_mean*N_data, axis=1)/np.nansum(N_data,axis=1)
OmA_stdv = np.nansum(OmA_stdv*N_data,axis=1)/np.nansum(N_data,axis=1)
Nobs_data = np.nansum(N_data, axis=1)

# Plotting

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
            
            lat_M36, lon_M36 = easev2_ind2latlon(np.arange(406), np.arange(964),'M36')
            lon_2d,lat_2d = np.meshgrid(lon_M36,lat_M36)
        else:
            grid_data, uy,ux = array2grid(tile_data, lat = tc['com_lat'], lon = tc['com_lon'])
            lon_2d,lat_2d = np.meshgrid(ux, uy)
            
        if 'normalized' in title_txt:
            title_txt = title_txt + '\n' + "avg=%.3f, avg(abs(nstdv-1))=%.3f" % (np.nanmean(grid_data), np.nanmean(np.abs(grid_data-1.)))+' '+units
        elif 'mean' in title_txt:
            title_txt = title_txt + '\n' + "avg=%.3f, avg(abs)=%.3f" % (np.nanmean(grid_data), np.nanmean(np.abs(grid_data)))+' '+units
        else:
            title_txt = title_txt + '\n' + "avg=%.2f" % (np.nanmean(grid_data)) +' '+units                
     
        if 'normalized' in title_txt:
            grid_data = np.log10(grid_data)
            crange = [-0.6, 0.45]
            
        mm, cs = plotMap(grid_data, ax =axes[i,j], lat=lat_2d, lon=lon_2d, cRange=crange, \
                    title=title_txt, cmap=colormap, bounding=[-60, 80, -180,180])            

plt.tight_layout()
# Save figure to file
fig.savefig(out_path+'Map_OmF_'+ expid +'_'+start_time.strftime('%Y%m')+'_'+\
                    end_time.strftime('%Y%m')+'.png')
#plt.show()
plt.close(fig)

