
import sys;       sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")
import os

import numpy as np

from datetime               import datetime, timedelta
from read_GEOSldas          import read_tilecoord, read_obs_param
from postproc_ObsFcstAna    import postproc_ObsFcstAna

def get_config():

    # =========================================================================================
    #
    # User-defined inputs for post-processing of ObsFcstAna output

    # Range of months to process:
    
    start_year  = 2015
    start_month =    4
    last_year   = 2016
    last_month  =    4

    # Sums or stats will be processed for exp_main:

    exp_main = {
        'expdir'        : '/discover/nobackup/projects/gmao/merra/iau/merra_land/SMAP_runs/SMAP_Nature_v11/',
        'expid'         : 'DAv8_SMOSSMAP',
        'exptag'        : 'DAMulti_SMAP', 
        'domain'        : 'SMAP_EASEv2_M36_GLOBAL',
        'da_t0'         : 3,                              # hour of first land analysis
        'da_dt'         : 10800,                          # ObsFcstAna file interval in seconds
        'species_list'  : [5,6,7,8],
        'obsparam_time' : "20150401_0000"                 # time stamp of obsparam file (YYYYMMDD_HHMM)
    }

    # Optional experiment(s) can be added for cross masking or extracting obs from a different experiment.
    #
    # All optional experiments and the main experiment must have identical tilecoords.
    # If the default "species" number/order do not match, set "species_list" accordingly to force a match.
    # Output will be cross-masked between all specified experiments.    

    # Forecasts and analyses are always from the main experiment.
    # Observations can be from experiment indicated by 'use_obs' set to True.
    # The mostly likely use case for this is that _scaled_ observations from a DA experiment
    # are used to compute OmF etc diagnostics for a corresponding open loop experiment.

    exp_sup1 = {
        'expdir'        : '/discover/nobackup/projects/gmao/merra/iau/merra_land/SMAP_runs/SMAP_Nature_v11/',
        'expid'         : 'DAv8_M36',
        'exptag'        : 'DASMAP_SMAP', 
        'domain'        : 'SMAP_EASEv2_M36_GLOBAL',
        'use_obs'     : True,
        'species_list'  : [1,2,3,4],
        'obsparam_time' : "20150401_0000"                 # time stamp of obsparam file (YYYYMMDD_HHMM)
    }

    # Convert experiments input to a list; first entry must be exp_main 

    exp_list = [exp_main, exp_sup1]
    
    # Top level directory to store monthly sum files; can use the experiment directory or a different path;
    # /Yyyyy/Mmm/ is added automatically for individual months

    out_path = '/discover/nobackup/qliu/SMAP_test/'
    
    sum_path = out_path+'/'+exp_main['expid']+'/output/'+exp_main['domain']+'/ana/ens_avg/'

    #
    # ===================== end of user-defined inputs =================================================

    # process time range info;  end_time is first of month after (end_year, end_month)

    if last_month==12 :
        end_month = 1
        end_year  = last_year + 1
    else :
        end_month = last_month + 1
        end_year  = last_year
        
    start_time = datetime( start_year, start_month, 1)
    end_time   = datetime( end_year,   end_month,   1)

    # Get tilecoord and obsparam information for each experiment
    
    for exp in exp_list:
        expdir   = exp['expdir']
        expid    = exp['expid']
        domain   = exp['domain']

        YYYY     = exp['obsparam_time'][0:4]
        MM       = exp['obsparam_time'][4:6]
        
        fop      = expdir+expid+'/output/'+domain+'/rc_out/Y'+YYYY+'/M'+MM+'/'+expid+'.ldas_obsparam.'+exp['obsparam_time']+'z.txt'
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
        tc  = read_tilecoord(ftc)

        exp.update({'tilecoord':tc, 'obsparam':obsparam})

    config ={
        'exp_list'   : exp_list,
        'start_time' : start_time,
        'end_time'   : end_time,
        'sum_path'   : sum_path,
        'out_path'   : out_path,
        }

    return config

# ====================== EOF =========================================================
