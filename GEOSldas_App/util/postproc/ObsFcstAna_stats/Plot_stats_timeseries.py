#!/usr/bin/env python3

"""
Sample script to compute spatially averaged monthly statistics
based on pre-saved monthly sums.

"""

import sys;       sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")
import os

import numpy             as np
import matplotlib.pyplot as plt

from datetime               import datetime, timedelta
from dateutil.relativedelta import relativedelta
from netCDF4                import Dataset, num2date

def get_mean_stats_from_sums(fnc4_sums):

    var_list = ['obs_obs', 'obs_fcst','obs_ana']

    mdata_sum = {}
    mdata2_sum = {}
    with Dataset(fnc4_sums,'r') as nc:
        mN_data = nc.variables['N_data'][:]
        moxf_sum = nc.variables['obsxfcst_sum'][:]
        moxa_sum = nc.variables['obsxana_sum'][:]
        moxf_sum[mN_data == 0] = np.nan
        moxa_sum[mN_data == 0] = np.nan
        for var in var_list:
            mdata_sum[var] = nc.variables[var+'_sum'][:]
            mdata2_sum[var] = nc.variables[var+'2_sum'][:]
            mdata_sum[var][mN_data == 0] = np.nan
            mdata2_sum[var][mN_data == 0] = np.nan

    # Make sure only aggregate tiles with valid values for all variables
    for var in var_list:
        mN_data = mN_data.astype(float)
        mN_data[np.isnan(mdata_sum[var])] = np.nan
        mN_data[mN_data == 0] = np.nan

    # cross mask before aggregating tile values   
    for var in var_list:
        mdata_sum[var][np.isnan(mN_data)] = np.nan
        mdata2_sum[var][np.isnan(mN_data)] = np.nan
        moxf_sum[np.isnan(mN_data)] = np.nan
        moxa_sum[np.isnan(mN_data)] = np.nan

    # Aggregate data of all tiles
    N_data = np.nansum(mN_data,axis=0)
    OxF_mean = np.nansum(moxf_sum,axis=0)/N_data
    OxA_mean = np.nansum(moxa_sum,axis=0)/N_data
    data_mean = {}
    data2_mean = {}
    data_var = {}
    for var in var_list:
        data_mean[var] = np.nansum(mdata_sum[var],axis=0)/N_data
        data2_mean[var] = np.nansum(mdata2_sum[var],axis=0)/N_data
        # var(x) = E[x2] - (E[x])^2
        data_var[var] = data2_mean[var] - data_mean[var]**2

    # Computer metrics of O-F, O-A, etc. based on above stats
    O_mean = data_mean['obs_obs']
    F_mean = data_mean['obs_fcst']
    A_mean = data_mean['obs_ana']
    O_var = data_var['obs_obs']
    F_var = data_var['obs_fcst']
    A_var = data_var['obs_ana']

    # mean(x-y) = E[x] - E[y]   
    OmF_mean = O_mean - F_mean
    OmA_mean = O_mean - A_mean

    # var(x-y) = var(x) + var(y) - 2cov(x,y)
    # cov(x,y) = E[xy] - E[x]E[y]
    OmF_stdv  = np.sqrt(O_var + F_var - 2 * (OxF_mean - O_mean*F_mean))
    OmA_stdv  = np.sqrt(O_var + A_var - 2 * (OxA_mean - O_mean*A_mean))

    # Combine metrics of individual species using weighted averaging
    Nobs_data = np.nansum(N_data)
    OmF_mean = np.nansum(OmF_mean*N_data)/Nobs_data
    OmF_stdv = np.nansum(OmF_stdv*N_data)/Nobs_data
    OmA_mean = np.nansum(OmA_mean*N_data)/Nobs_data
    OmA_stdv = np.nansum(OmA_stdv*N_data)/Nobs_data

    return OmF_mean, OmF_stdv, OmA_mean, OmA_stdv, Nobs_data
    
if __name__ == "__main__":

    import calendar
    import pickle

    # Define time range for processing
    start_time = datetime(2015,4,1)
    end_time = datetime(2016,4,1)

    expid = 'Vv7032'
    mo_sums_path = '/discover/nobackup/qliu/SMAP_diag/SPL4SM_Vv7032/'+ \
              'output/SMAP_EASEv2_M09_GLOBAL/ana/ens_avg/'

    date_time = start_time
    OmF_mean =[]
    OmF_stdv =[]
    OmA_mean =[]
    OmA_stdv=[]
    Ndata =[]
    date_vec =[]

    while date_time < end_time:
        print('compute monthly spatial mean stats for '+date_time.strftime('%Y%m'))
        fname_sums = mo_sums_path+'/Y'+date_time.strftime('%Y')+ \
            '/M'+date_time.strftime('%m')+'/'+expid+ \
            '.ens_avg.ldas_ObsFcstAna.'+date_time.strftime('%Y%m')+'_sums.nc4'

        if os.path.isfile(fname_sums):
            OmFm,OmFs,OmAm,OmAs,Nobsm = get_mean_stats_from_sums(fname_sums)
        else:
            print(f"ERROR: File '{fname_sums}' not found.")
            sys.exit(1)
        
        OmF_mean.append(OmFm)
        OmF_stdv.append(OmFs)
        OmA_mean.append(OmAm)
        OmA_stdv.append(OmAs)
        Ndata.append(Nobsm)

        date_vec.append(date_time.strftime('%Y%m'))
        date_time = date_time + relativedelta(months=1)

    stats_dict = {"OmF_mean":OmF_mean,
                         "OmF_stdv":OmF_stdv,
                          "OmA_mean":OmA_mean,
                          "OmA_stdv":OmA_stdv,
                          "Ndata": Ndata}
    
    plot_var = 'OmF_mean'

    fig, ax = plt.subplots(figsize=(10,4))
    bars = ax.bar(date_vec, stats_dict[plot_var])

    plt.grid(True, linestyle='--', alpha=0.5)

    plt.xticks(ticks=date_vec[::2], labels=date_vec[::2])
    plt.title(expid+ 'monthly '+plot_var)
    plt.xlim(-1, len(date_vec)+1)
    plt.ylim(-.7, 1.5)

    plt.tight_layout()
    plt.show()
    plt.savefig('Bars_'+plot_var+expid+start_time.strftime('%Y%m')+'_'+\
                    end_time.strftime('%Y%m')+'.png')
    plt.close()

# ====================== EOF =========================================================
