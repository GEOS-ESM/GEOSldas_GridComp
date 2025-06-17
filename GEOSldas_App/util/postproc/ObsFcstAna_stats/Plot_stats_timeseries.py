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

from datetime               import datetime, timedelta
from dateutil.relativedelta import relativedelta
    
def Plot_monthly_OmF_bars(postproc_obj, fig_path='./'):
    import pickle
    exptag     = postproc_obj.exptag
    start_time = postproc_obj.start_time
    end_time   = postproc_obj.end_time

    stats_file = fig_path + 'spatial_stats_'+exptag+'_'+start_time.strftime('%Y%m')+ \
        '_'+(end_time+timedelta(days=-1)).strftime('%Y%m')+'.pkl'
    
    if os.path.isfile(stats_file):
        
        with open(stats_file,'rb') as file:
            stats = pickle.load(file)
        date_vec = stats['date_vec']
        
    else:
        
        # Initialize monthly metrics as list
        OmF_mean =[]
        OmF_stdv =[]
        OmA_mean =[]
        OmA_stdv =[]
        Ndata    =[]
        date_vec =[]

        # Time loop 
        current_time = start_time
        while current_time < end_time:
            print('compute monthly spatial mean stats for '+ current_time.strftime('%Y%m'))

            stats_mo = postproc_obj.calc_spatial_stats_from_sums(current_time)

            OmFm = stats_mo['OmF_mean']
            OmFs = stats_mo['OmF_stdv']
            OmAm = stats_mo['OmA_mean']
            OmAs = stats_mo['OmA_stdv']
            Nobsm = stats_mo['N_data']
                  
            Ndata.append(Nobsm)

            OmF_mean.append(OmFm)
            OmF_stdv.append(OmFs)
            OmA_mean.append(OmAm)
            OmA_stdv.append(OmAs)

            date_vec.append(current_time.strftime('%Y%m'))
            current_time = current_time + relativedelta(months=1)

        # Store stats in a dictionary for easier saving and referencing  
        stats = {"OmF_mean":np.array(OmF_mean), "OmF_stdv":np.array(OmF_stdv),
                 "OmA_mean":np.array(OmA_mean), "OmA_stdv":np.array(OmA_stdv),
                 "Ndata":np.array(Ndata), "date_vec":date_vec}
        
        if stats_file is not None:
            with open(stats_file,'wb') as file:
                pickle.dump(stats,file)

    # Compute Nobs-weighted avg of each metric across all species.
    # Typically used for SMAP Tb_h/h from asc and desc overpasses,
    # or ASCAT soil moisture from Metop-A/B/C.
    # DOES NOT MAKE SENSE IF, SAY, SPECIES HAVE DIFFERENT UNITS!
    stats_plot = {}
    Ndata = np.nansum(stats['Ndata'], axis=1)
    stats_plot['Ndata']     = Ndata
    stats_plot['OmF_mean']  = np.nansum(stats['OmF_mean']*stats['Ndata'], axis=1)/Ndata
    stats_plot['OmF_stdv']  = np.nansum(stats['OmF_stdv']*stats['Ndata'], axis=1)/Ndata
    stats_plot['OmA_mean']  = np.nansum(stats['OmA_mean']*stats['Ndata'], axis=1)/Ndata
    stats_plot['OmA_stdv']  = np.nansum(stats['OmA_stdv']*stats['Ndata'], axis=1)/Ndata
    
    plot_var = 'OmF_mean'

    fig, ax = plt.subplots(figsize=(10,4))
    bars = ax.bar(date_vec, stats_plot[plot_var])

    plt.grid(True, linestyle='--', alpha=0.5)

    plt.xticks(ticks=date_vec[::2], labels=date_vec[::2])
    plt.title(exptag+ ' monthly '+plot_var)
    plt.xlim(-1, len(date_vec)+1)
    plt.ylim(-.1, 2.)

    plt.tight_layout()
    #plt.show()
    plt.savefig(fig_path+'Bars_'+plot_var+exptag+date_vec[0]+'_'+\
                    date_vec[-1]+'.png')
    plt.close()

if __name__ == "__main__":

    from postproc_ObsFcstAna    import postproc_ObsFcstAna
    from user_config            import get_config

    config     = get_config()                 # edit user-defined inputs in user_config.py
    
    exp_list   = config['exp_list']
    start_time = config['start_time']
    end_time   = config['end_time']
    sum_path   = config['sum_path']
    out_path   = config['out_path']
     
    # ------------------------------------------------------------------------------------
    #
    # Initialize the postprocessing object
    postproc = postproc_ObsFcstAna(exp_list, start_time, end_time, sum_path=sum_path)

    # Compute spatial stats and plot monthly O-F stats bars
    Plot_monthly_OmF_bars(postproc, fig_path=out_path)


# ====================== EOF =========================================================
