#!/usr/bin/env python3

"""
Sample script for plotting monthly time series of spatially averaged data assimilation diagnostics.
Computes Nobs-weighted avg of each metric across all species.
Requires saved files with monthly sums (see Get_ObsFcstAna_stat.py).
"""

import sys;       sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")
import os

import numpy             as np
import matplotlib.pyplot as plt
import pickle

from datetime               import datetime, timedelta
from dateutil.relativedelta import relativedelta

from postproc_ObsFcstAna    import postproc_ObsFcstAna
from user_config            import get_config
    
def Main_OmF_timeseries():

    config     = get_config()                 # edit user-defined inputs in user_config.py
    
    exp_list   = config['exp_list']
    start_time = config['start_time']
    end_time   = config['end_time']
    sum_path   = config['sum_path']
    out_path   = config['out_path']
     
    # ------------------------------------------------------------------------------------
    #
    # Initialize the postprocessing object
    postproc_obj = postproc_ObsFcstAna(exp_list, start_time, end_time, sum_path=sum_path)
    
    exptag     = postproc_obj.exptag
    start_time = postproc_obj.start_time
    end_time   = postproc_obj.end_time

    stats_file = out_path + 'spatial_stats_'+exptag+'_'+start_time.strftime('%Y%m')+ \
        '_'+(end_time+timedelta(days=-1)).strftime('%Y%m')+'.pkl'
    
    if os.path.isfile(stats_file):
        
        with open(stats_file,'rb') as file:
            print(f'reading from {stats_file}')
            stats = pickle.load(file)
        
    else:

        stats = postproc_obj.calc_spatial_stats_from_sums()
        
        if stats_file is not None:
            with open(stats_file,'wb') as file:
                print(f'saveing stats to {stats_file}')
                pickle.dump(stats,file)

    date_vec = stats['date_vec']

    # Compute Nobs-weighted avg of each metric across all species.
    # Typically used for SMAP Tb_h/h from asc and desc overpasses,
    # or ASCAT soil moisture from Metop-A/B/C.
    # DOES NOT MAKE SENSE IF, SAY, SPECIES HAVE DIFFERENT UNITS!
    stats_plot = {}
    N_data = np.nansum(stats['N_data'], axis=1)
    stats_plot['N_data']     = N_data
    stats_plot['OmF_mean']  = np.nansum(stats['OmF_mean']*stats['N_data'], axis=1)/N_data
    stats_plot['OmF_stdv']  = np.nansum(stats['OmF_stdv']*stats['N_data'], axis=1)/N_data
    stats_plot['OmA_mean']  = np.nansum(stats['OmA_mean']*stats['N_data'], axis=1)/N_data
    stats_plot['OmA_stdv']  = np.nansum(stats['OmA_stdv']*stats['N_data'], axis=1)/N_data
    
    plot_var = 'OmF_mean'

    fig, ax = plt.subplots(figsize=(10,4))
    bars = ax.bar(date_vec, stats_plot[plot_var])

    plt.grid(True, linestyle='--', alpha=0.5)

    plt.xticks(ticks=date_vec[::2], labels=date_vec[::2])
    plt.title(exptag+ 'monthly '+plot_var)
    plt.xlim(-1, len(date_vec)+1)
    plt.ylim(-.1, 2.)

    plt.tight_layout()
    #plt.show()
    plt.savefig(out_path+'Bars_'+plot_var+exptag+date_vec[0]+'_'+\
                    date_vec[-1]+'.png')
    plt.close()

if __name__ == "__main__":
    
    Main_OmF_timeseries()


# ====================== EOF =========================================================
