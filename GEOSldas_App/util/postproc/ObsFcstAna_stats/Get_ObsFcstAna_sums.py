#!/usr/bin/env python3

"""
Sample script for post-processing GEOSldas ObsFcstAna output into data assimilation diagnostics.
First, compute and store monthly sums and sums of squares and cross-products of raw ObsFcstAna output.  
Data assimilation diagnostics ("stats") such as the mean and std-dev of the observation-minus-forecast 
residuals can then be diagnosed quickly from these intermediate "sums" files. 
Sample script optionally computes and plots:
- Maps of long-term data assimilation diagnostics (see also Plot_stats_maps.py).
- Monthly time series of spatially averaged data assimilation diagnostics (see also Plot_stats_timeseries.py).

Usage:
    1. Edit "user_config.py" with experiments information.
    2. Run this script as follows (on Discover):
    $ module load python/GEOSpyD
    $ ./Get_ObsFcstAna_sums.py
    
    # Background execution:
    $ nohup ./Get_ObsFcstAna_sums.py > out.log &

Authors: Q. Liu, R. Reichle, A. Fox
Last Modified: May 2025
"""

import sys;       sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")
import os

import numpy as np

from datetime               import timedelta
from postproc_ObsFcstAna    import postproc_ObsFcstAna
from user_config            import get_config              # user-defined config info

# ---
#
# If the script is run in the background, uncomment the following lines to see the redirected
#   standard output in the out.log file immediately.  When the lines are commented out, the redirected
#   standard output will not appear in the out.log file until the job has completed.
# If the script is run in the foreground, the lines can be commented out to monitor the output log.
#
#import io
#sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
#sys.stderr = io.TextIOWrapper(open(sys.stderr.fileno(), 'wb', 0), write_through=True)
#
# ---

def main():
     
    config     = get_config()                 # edit user-defined inputs in user_config.py
    
    exp_list   = config['exp_list']
    start_time = config['start_time']
    end_time   = config['end_time']
    sum_path   = config['sum_path']
    out_path   = config['out_path']
     
    # ------------------------------------------------------------------------------------
    # Postprocess raw ObsFcstAna output data into monthly sums 

    # Initialize the postprocessing object
    postproc = postproc_ObsFcstAna(exp_list, start_time, end_time, sum_path=sum_path)

    # Compute and save monthly sums
    postproc.save_monthly_sums()

if __name__ == '__main__':
    main()

# ====================== EOF =========================================================
