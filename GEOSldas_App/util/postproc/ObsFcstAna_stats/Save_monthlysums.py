#!/usr/bin/env python3

"""
Sample script for computing and storing monthly sums and sums of squares and cross-products
from GEOSldas ObsFcstAna output.  Data assimilation diagnostics such as the mean and
std-dev of the observation-minus-forecast residuals can be diagnosed quickly from these
intermediate "sums" files. 

Usage on Discover:
    First, edit "user_config.py" with proper experiments inforamtion.
    Then, run this script as follows
    $ module load python/GEOSpyD
    $ ./Save_monthlysums.py
    
    # Background execution:
    $ nohup ./Save_monthlysums.py > out.log &

Authors: Q. Liu, R. Reichle, A. Fox
Last Modified: May 2025
"""

import sys;       sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")
import os

import numpy as np

from datetime               import datetime, timedelta
from read_GEOSldas          import read_tilecoord, read_obs_param
from postproc_ObsFcstAna    import postproc_ObsFcstAna
from user_config            import get_config
# If the script is run in the background, uncomment the following lines to see the redirected
#   standard output in the out.log file immediately.  When the lines are commented out, the redirected
#   standard output will not appear in the out.log file until the job has completed.
# If the script is run in the foreground, the lines must be commented out.
#
#import io
#sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
#sys.stderr = io.TextIOWrapper(open(sys.stderr.fileno(), 'wb', 0), write_through=True)


def main():
     
    config = get_config()
    
    exp_list = config['exp_list']
    start_time = config['start_time']
    end_time = config['end_time']
    obs_from = config['obs_from']
    sum_path = config['sum_path']
    out_path = config['out_path']
     
    # ------------------------------------------------------------------------------------
    # First
    # Postprocess raw ObsFcstAna output data into monthly sums 

    # Initialize the postprocessing object
    postproc = postproc_ObsFcstAna(exp_list, start_time, end_time, obs_from=obs_from, sum_path=sum_path)

    # Compute and save monthly sums
    postproc.save_monthly_sums()

    # --------------------------------------------------------------------------------------
    # Next, there is an option to compute long-term temporal/spatial statistics and create plots.
    # The plotting scripts can also run standalone using the individual Plot_xx.py scripts,
    # as long as the monthly sum files are available.
    
    plot_maps = False
    plot_timeseries = False

    if plot_maps:
        # Compute long-term temporal stats and plot maps
        if len(exp_list) >1 :
            stats_file  = out_path + 'tmp_stats_'+exp_list[0]['exptag']+'_obsfrom_'+ \
                          exp_list[obs_from]['exptag']+'_'+start_time.strftime('%Y%m%d')+'_'+ \
                          (end_time+timedelta(days=-1)).strftime('%Y%m%d')+'.nc4'
        else:
            stats_file  = out_path + 'tmp_stats_'+exp_list[0]['exptag']+'_'+ start_time.strftime('%Y%m%d')+'_'+ \
                          (end_time+timedelta(days=-1)).strftime('%Y%m%d')+'.nc4'

        # temporal_stats is a dictionary that contains all mean/variances fields for computing long-term O-F/O-A stats
        # each field has the dimension [N_tile, N_species]

        temporal_stats = postproc.calc_temporal_stats_from_sums(write_to_nc=True, fout_stats=stats_file)

        # Example to plot some O-F maps
        from Plot_stats_maps import plot_OmF_maps
        plot_OmF_maps(postproc, temporal_stats, fig_path=out_path )


    if plot_timeseries:
        # Example to compute spatial stats and plot monthly O-F stats bars
        from Plot_stats_timeseries import Plot_monthly_OmF_bars
        Plot_monthly_OmF_bars(postproc, fig_path=out_path)

if __name__ == '__main__':
    main()

# ====================== EOF =========================================================
