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
    
def Plot_monthly_OmF_bars(postproc_obj, fig_path='./'):
    import pickle
    expid = postproc_obj.outid
    start_month = postproc_obj.start_time
    end_month = postproc_obj.end_time

    stats_file = fig_path + 'spa_stats_'+expid+'_'+start_month.strftime('%Y%m')+ \
            '_'+(end_month+timedelta(days=-1)).strftime('%Y%m')+'.pkl'
    
    if os.path.isfile(stats_file):
        
        with open(stats_file,'rb') as file:
            stats_dict = pickle.load(file)
        date_vec = stats_dict['date_vec']
        
    else:
        
        # Initialize monthly metrics as list
        OmF_mean =[]
        OmF_stdv =[]
        OmA_mean =[]
        OmA_stdv=[]
        Ndata =[]
        date_vec =[]

        # Time loop 
        current_month = start_month
        while current_month < end_month:
            print('compute monthly spatial mean stats for '+ current_month.strftime('%Y%m'))

            OmFm,OmFs,OmAm,OmAs,Nobsm = postproc_obj.calc_spatial_stats_from_sums(current_month)

            # Average individual species into a single value
            Nobsm = np.nansum(Nobsm)
            OmFm = np.nansum(OmFm*Nobsm)/Nobsm
            OmFs = np.nansum(OmFs*Nobsm)/Nobsm
            OmAm = np.nansum(OmAm*Nobsm)/Nobsm
            OmAs = np.nansum(OmAs*Nobsm)/Nobsm
      
            OmF_mean.append(OmFm)
            OmF_stdv.append(OmFs)
            OmA_mean.append(OmAm)
            OmA_stdv.append(OmAs)
            Ndata.append(Nobsm)

            date_vec.append(current_month.strftime('%Y%m'))
            current_month = current_month + relativedelta(months=1)

        # Store stats in a dictionary for easier saving and referencing  
        stats_dict = {"OmF_mean":OmF_mean, "OmF_stdv":OmF_stdv,
                             "OmA_mean":OmA_mean, "OmA_stdv":OmA_stdv,
                             "Ndata": Ndata, "date_vec":date_vec}
        
        if stats_file is not None:
            with open(stats_file,'wb') as file:
                pickle.dump(stats_dict,file)
    
    plot_var = 'OmF_mean'

    fig, ax = plt.subplots(figsize=(10,4))
    bars = ax.bar(date_vec, stats_dict[plot_var])

    plt.grid(True, linestyle='--', alpha=0.5)

    plt.xticks(ticks=date_vec[::2], labels=date_vec[::2])
    plt.title(expid+ 'monthly '+plot_var)
    plt.xlim(-1, len(date_vec)+1)
    plt.ylim(-.1, 2.)

    plt.tight_layout()
    #plt.show()
    plt.savefig(fig_path+'Bars_'+plot_var+expid+date_vec[0]+'_'+\
                    date_vec[-1]+'.png')
    plt.close()

if __name__ == "__main__":
    from postproc_ObsFcstAna    import postproc_ObsFcstAna
    from user_config            import get_config
    config = get_config()
    
    exp_list = config['exp_list']
    start_time = config['start_time']
    end_time = config['end_time']
    sum_path = config['sum_path']
    out_path = config['out_path']
     
    # ------------------------------------------------------------------------------------
    # First
    # Postprocess raw ObsFcstAna output data into monthly sums 

    # Initialize the postprocessing object
    postproc = postproc_ObsFcstAna(exp_list, start_time, end_time, sum_path=sum_path)

    Plot_monthly_OmF_bars(postproc, fig_path=out_path)


# ====================== EOF =========================================================
