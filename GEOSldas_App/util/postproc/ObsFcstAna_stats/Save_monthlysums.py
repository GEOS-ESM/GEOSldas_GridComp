#!/usr/bin/env python3

"""
Sample script for processing data assimilation (DA) diagnostics from 
GEOSldas ObsFcstAna output into monthly sums. Long-term statistics can be computed
efficiently later using the monthly data. 

Processing workflow:
1. User provides experiment information and calls postproc_ObsFcstAna to:
   - Calculate/save monthly sums of Obs/Fcst/Ana and squared/products
   - Calculate mean/stdv of Obs/Fcst/Ana
2. User computes desired DA diagnostic statistics and creates plots

Usage on Discover:
    $ module load python/GEOSpyD
    $ ./Main_example.py
    
    # Background execution:
    $ nohup ./Main_example.py > out.log &

Author: Q. Liu, R. Reichle, A. Fox
Last Modified: May., 2025
"""

import sys;       sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")
import os

import numpy             as np
from datetime               import datetime, timedelta

from read_GEOSldas          import read_tilecoord, read_obs_param

from postproc_ObsFcstAna    import postproc_ObsFcstAna

# Uncomment if to run the script in the background to see the standard output while running 
# import io
#sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
#sys.stderr = io.TextIOWrapper(open(sys.stderr.fileno(), 'wb', 0), write_through=True)

# User provided time range for processing (Year, Month, Day)
start_time = datetime(2016,4,1)
end_time   = datetime(2016,5,1)

# -------------------------------- Experiments Information -----------------------------------------
# Supports single experiment or multiple experiments.
# All experiments must have identical tilecoords and number/order of observation species.
# If the default "species" number/order do not match, need to set the *optional*
#   "select_species" key to get a match, i.e. same species sequences.
# This capability is required to enable calculating OmF/OmA statistics for one experiment
#   using observations from another experiment. See note below.

# User provide at lease one experiment 
exp_main = { 'expdir'       : '/discover/nobackup/projects/gmao/merra/iau/merra_land/SMAP_runs/SMAP_Nature_v11/',
          'expid'        : 'DAv8_SMOSSMAP',
          'exptag'       : 'DAMulti_SMAP', 
          'domain'       : 'SMAP_EASEv2_M36_GLOBAL',
          'da_t0'        : 3,       # first hour of the month 
          'da_dt'        : 10800,   # ObsFcstAna file interval in seconds
          'species_list' : [5,6,7,8] }

# Optional experiment can be added for cross masking or extracting obs from a different experiment
exp_sup1 = { 'expdir'       : '/discover/nobackup/projects/gmao/merra/iau/merra_land/SMAP_runs/SMAP_Nature_v11/',
          'expid'        : 'DAv8_M36',
          'exptag'       : 'DASMAP_SMAP', 
          'domain'       : 'SMAP_EASEv2_M36_GLOBAL',
          'species_list' : [1,2,3,4]  }

# User provided directory to store monthly sum files
# can use the main experiment's ana/ or a different path
monthly_sums_path = '/discover/nobackup/qliu/SMAP_diag/' +exp_main['expid']+ \
              '/output/'+exp_main['domain']+'/ana/ens_avg/'

# forecasts and analyses are always from the main experiment.
# observations can be from experiment indicated by 'obs_from' index.
# The mostly likely use case for this is that _scaled_ observations from a DA experiment
# are used to compute OmF etc diagnostics for a corresponding open loop experiment.

exp_list = [exp_main, exp_sup1]
obs_from = 1                            # obs is from "exp_sup1" (0-based indexing)
if obs_from >= len(exp_list):
    print('Invalid "obs_from" value')
    sys.exit()

# Read tilecoord and obs_param information of each experiment
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

#  =========================================================================
#  Postprocess raw ObsFcstAna output data into monthly sums for simpler and faster postprocessing;
#  computes mean, variance from monthly sums that can be used to compute DA diagnostics directly

# Initialize the postprocessing object
postproc = postproc_ObsFcstAna(exp_list, start_time, end_time, obs_from=obs_from)
# Compute and save monthly sums 
postproc.save_monthly_sum(monthly_sums_path)

print('Processing completed!')   
# ====================== EOF =========================================================
