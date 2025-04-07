#!/usr/bin/env python3

"""
GEOSldas DA Diagnostics Map Generator

This script computes and plots GEOSldas Data Assimilation (DA) diagnostics maps
based on the FcstObsAna output. It processes monthly data and generates final
metrics from these monthly statistics.

The script performs the following main tasks:
1. Loads and processes monthly FcstObsAna output data
2. Computes various DA diagnostic metrics
3. Generates maps and plots of the computed metrics
4. Saves temporary monthly statistics and final aggregated metrics

Usage:
    On NASA's Discover:
    $ module load python/GEOSpyD
    $ ./main_ObsFcstAna.py
     or to run in the background,
    $ nohup ./main_ObsFcstAna.py > out.log &

Requirements:
    - Python 3.x
    - Modules: numpy, matplotlib, netCDF4 (included in GEOSpyD)

Author: Q. Liu
Last Modified: Apr 1, 2025
"""

import numpy as np
import os
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from helper.read_GEOSldas import read_tilecoord, read_obs_param
from helper.util import make_folder, array2grid
from helper.plot import plotMap
from helper.smapeasev2 import smapeasev2_ind2latlon
from helper.compute_monthly_stats import compute_monthly_stats
from helper.write_nc4 import write_sums_nc4

import warnings; warnings.filterwarnings("ignore")
import sys 
import io

sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
sys.stderr = io.TextIOWrapper(open(sys.stderr.fileno(), 'wb', 0), write_through=True)

expdir = '/gpfsm/dnb05/projects/p51/SMAP_Nature/'
expid = 'SPL4SM_Vv8010'
domain = 'SMAP_EASEv2_M09_GLOBAL'

start_time = datetime(2015,4,1)
end_time = datetime(2021,4,1)

# Define a minimum threshold for the temporal data points to ensure statistical reliability
# of the computed metrics. 
Nmin = 20

# Base directory for storing monthly files
# This can be the same as the experiment directory (expdir) or a different location
out_path_mo = '/discover/nobackup/qliu/SMAP_diag/' +expid+'/output/'+domain+'/ana/ens_avg/'

# Directory for diagnostic plots
out_path = '/discover/nobackup/qliu/SMAP_diag/'
make_folder(out_path)

# Variable list for computing sum and sum of squared
var_list = ['obs_obs', 'obs_obsvar','obs_fcst','obs_fcstvar','obs_ana','obs_anavar']

# Read tilecoord and obsparam for tile and obs species information
ftc = expdir+expid+'/output/'+domain+'/rc_out/'+expid+'.ldas_tilecoord.bin'
tc = read_tilecoord(ftc)
n_tile = tc['N_tile']

fop = expdir+expid+'/output/'+domain+'/rc_out/Y2015/M04/'+expid+'.ldas_obsparam.20150401_0000z.txt'
obs_param = read_obs_param(fop)
n_spec = len(obs_param)

# Initialize statistical metrics 
data_sum = {}
data2_sum = {}
N_data = np.zeros((n_tile, n_spec))
oxf_sum = np.zeros((n_tile, n_spec))
oxa_sum = np.zeros((n_tile, n_spec))
fxa_sum = np.zeros((n_tile, n_spec))

for var in var_list:
    data_sum[var] = np.zeros((n_tile, n_spec))
    data2_sum[var] = np.zeros((n_tile, n_spec))

# Time loop: processing data at monthly time step
date_time = start_time
while date_time < end_time:
    # File to store monthly statistics    
    fout_path = out_path_mo + '/Y'+ date_time.strftime('%Y') + '/M' + date_time.strftime('%m') + '/'
    make_folder(fout_path)
    
    fout = fout_path + expid+'.ens_avg.ldas_ObsFcstAna.' + date_time.strftime('%Y%m') +'_stats.nc4'

    # Read monthly data if file exists, otherwise compute monthly statistics first   
    if os.path.isfile(fout):
        print('read sums from  monthly file: '+fout)
        mdata_sum = {}
        mdata2_sum = {}
        with Dataset(fout,'r') as nc:
            mN_data = nc.variables['N_data'][:]
            moxf_sum = nc.variables['obsxfcst_sum'][:]
            moxa_sum = nc.variables['obsxana_sum'][:]
            mfxa_sum = nc.variables['fcstxana_sum'][:]
            for var in var_list:
                mdata_sum[var] = nc.variables[var+'_sum'][:]
                mdata2_sum[var] = nc.variables[var+'2_sum'][:]
    else:
        print('compute monthly sums for '+date_time.strftime('%Y%m'))
        mN_data, mdata_sum, mdata2_sum, moxf_sum, moxa_sum, mfxa_sum = \
                 compute_monthly_stats(expdir,expid,domain,date_time,tc,obs_param,var_list)
        print('save to monthly file: '+fout)
        write_sums_nc4(fout, mN_data,mdata_sum, mdata2_sum, moxf_sum, moxa_sum, mfxa_sum, obs_param)

    # Aggregate monthly data
    N_data += mN_data
    oxf_sum += moxf_sum
    oxa_sum += moxa_sum
    fxa_sum += mfxa_sum
   
    for var in var_list:
        data_sum[var] += mdata_sum[var] 
        data2_sum[var] += mdata2_sum[var]  
        
    date_time =date_time + relativedelta(months=1)

# Compute the final statistics
# This section calculate the final statistical metrics based on the accumulated data.
data_mean ={}
data2_mean = {}
data_var = {}

# First, compute the metrics of individual variables  
for var in var_list:
    data_sum[var][N_data == 0] = np.nan
    data2_sum[var][N_data == 0] = np.nan
    
    data_mean[var]  = data_sum[var] / N_data
    data2_mean[var] = data2_sum[var] /N_data
    # var(x) = E[x2] - (E[x])^2
    data_var[var] = data2_mean[var] - data_mean[var]**2
    
oxf_sum[N_data == 0] = np.nan
oxa_sum[N_data == 0] = np.nan
fxa_sum[N_data == 0] = np.nan
# E[xy]
oxf_mean = oxf_sum / N_data
oxa_mean = oxa_sum / N_data
fxa_mean = fxa_sum / N_data

# Then computer metrics of O-F, O-A, etc. based on above computed 
# mean(x-y) = E[x] - E[y]   
OmF_mean = data_mean['obs_obs'] - data_mean['obs_fcst']
OmA_mean = data_mean['obs_obs'] - data_mean['obs_ana']
# var(x-y) = var(x) + var(y) - 2cov(x,y)
# cov(x,y) = E[xy] - E[x]E[y]
OmF_stdv  = np.sqrt(data_var['obs_obs'] + data_var['obs_fcst'] - \
                       2 * (oxf_mean - data_mean['obs_obs']*data_mean['obs_fcst']))
                    
OmA_stdv  = np.sqrt(data_var['obs_obs'] + data_var['obs_ana'] - \
                       2 * (oxa_mean - data_mean['obs_obs']*data_mean['obs_ana']))

OmF_norm_mean = OmF_mean / np.sqrt(data_mean['obs_obsvar'] + data_mean['obs_fcstvar']) 
OmF_norm_stdv = np.sqrt(OmF_stdv**2 / (data_mean['obs_obsvar'] + data_mean['obs_fcstvar']) )
    
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
            if i==0 and j==0:
                grid_data = np.sum(reshaped,axis=(1, 3)) 
            else:
                grid_data = np.nanmean(reshaped,axis=(1, 3))
                
            lat_M36, lon_M36 = smapeasev2_ind2latlon(np.arange(406), np.arange(964),'M36')
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
fig.savefig(out_path+'Map_OmF_'+expid+'_'+start_time.strftime('%Y%m')+'_'+\
                    end_time.strftime('%Y%m')+'.png')
#plt.show()
plt.close(fig)

