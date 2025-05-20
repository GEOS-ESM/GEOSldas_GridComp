#!/usr/bin/env python3

"""
Sample script for computing and storing monthly sums and sums of squares and cross-products
from GEOSldas ObsFcstAna output.  Data assimilation diagnostics such as the mean and
std-dev of the observation-minus-forecast residuals can be diagnosed quickly from these
intermediate "sums" files.

Usage on Discover:
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

# If the script is run in the background, uncomment the following lines to see the redirected
#   standard output in the out.log file immediately.  When the lines are commented out, the redirected
#   standard output will not appear in the out.log file until the job has completed.
# If the script is run in the foreground, the lines must be commented out.
#
#import io
#sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
#sys.stderr = io.TextIOWrapper(open(sys.stderr.fileno(), 'wb', 0), write_through=True)

# --------------------------------------------------------------------------------------------------
#
# User-defined inputs

# Range of months to process:

start_year  = 2016
start_month =    4
last_year   = 2016
last_month  =    4

# Experiments to process:

# Supports single experiment or multiple experiments.
# All experiments must have identical tilecoords and number/order of observation species.
# If the default "species" number/order do not match, need to set the *optional*
#   "select_species" key to get a match, i.e. same species sequences.
# This capability is required to enable calculating OmF/OmA statistics for one experiment
#   using observations from another experiment. See note below.

# User must provide at least one experiment 

exp_main = {
    'expdir'        : '/discover/nobackup/projects/gmao/merra/iau/merra_land/SMAP_runs/SMAP_Nature_v11/',
    'expid'         : 'DAv8_SMOSSMAP',
    'exptag'        : 'DAMulti_SMAP', 
    'domain'        : 'SMAP_EASEv2_M36_GLOBAL',
    'da_t0'         : 3,                              # hour of first land analysis
    'da_dt'         : 10800,                          # ObsFcstAna file interval in seconds
    'species_list'  : [5,6,7,8],
    'obsparam_time' : "20150401_0000"                 # YYYYMMDD_HHMM
}

# Optional experiment can be added for cross masking or extracting obs from a different experiment

exp_sup1 = {
    'expdir'        : '/discover/nobackup/projects/gmao/merra/iau/merra_land/SMAP_runs/SMAP_Nature_v11/',
    'expid'         : 'DAv8_M36',
    'exptag'        : 'DASMAP_SMAP', 
    'domain'        : 'SMAP_EASEv2_M36_GLOBAL',
    'species_list'  : [1,2,3,4],
    'obsparam_time' : "20150401_0000"                 # YYYYMMDD_HHMM    
}

# User provided directory to store monthly sum files;  can use the main experiment's ana/ or a different path

out_path = '/discover/nobackup/qliu/SMAP_diag/'

# forecasts and analyses are always from the main experiment.
# observations can be from experiment indicated by 'obs_from' index.
# The mostly likely use case for this is that _scaled_ observations from a DA experiment
# are used to compute OmF etc diagnostics for a corresponding open loop experiment.

obs_from = 1              # take obs from "exp_sup1" (0-based indexing)

# ---------------------------------------------------------------------------------------------------------

# process time range info;  end_time is first of month after (end_year, end_month)

if last_month==12 :
    end_month = 1
    end_year  = last_year + 1
else :
    end_month = last_month + 1
    end_year  = last_year
    
start_time = datetime( start_year, start_month, 1)
end_time   = datetime( end_year,   end_month,   1)

# assemble output path

monthly_sums_path = out_path + '/' + exp_main['expid']+ '/output/' + exp_main['domain'] + '/ana/ens_avg/'

# process experiment info

exp_list = [exp_main, exp_sup1]

if obs_from >= len(exp_list):
    print('Invalid "obs_from" value')
    sys.exit()

# Read tilecoord and obs_param information of each experiment

for exp in exp_list:
    expdir   = exp['expdir']
    expid    = exp['expid']
    domain   = exp['domain']

    YYYY     = exp['obsparam_time'][0:4]
    MM       = exp['obsparam_time'][4:6]
    
    fop      = expdir+expid+'/output/'+domain+'/rc_out/Y'+YYYY+'/M'+MM+'/'+expid+'.ldas_obsparam.'+exp[obsparam_time]+'z.txt'
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

# Postprocess raw ObsFcstAna output data into monthly sums 

# Initialize the postprocessing object
postproc = postproc_ObsFcstAna(exp_list, start_time, end_time, obs_from=obs_from)

# Compute and save monthly sums
postproc.save_monthly_sum(monthly_sums_path)

print('Save_monthlysums: Processing completed!')   

# ====================== EOF =========================================================
