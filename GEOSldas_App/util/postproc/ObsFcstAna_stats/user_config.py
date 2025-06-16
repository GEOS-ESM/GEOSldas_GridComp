
import sys;       sys.path.append('../../shared/python/')
import warnings;  warnings.filterwarnings("ignore")
import os

import numpy as np

from datetime               import datetime, timedelta
from read_GEOSldas          import read_tilecoord, read_tilegrids, read_obs_param
from postproc_ObsFcstAna    import postproc_ObsFcstAna

def get_config():

    # =========================================================================================
    #
    # User-defined inputs for post-processing of ObsFcstAna output

    # Range of months to process:
    
    start_year  = 2015
    start_month =    4
    last_year   = 2016
    last_month  =    3

    # Sums or stats will be processed for exp_main:

    exp_main = {
        'expdir'        : '/discover/nobackup/projects/gmao/merra/iau/merra_land/SMAP_runs/SMAP_Nature_v11/',
        'expid'         : 'DAv8_SMOSSMAP',             # GEOSldas exp ID of simulation
        'exptag'        : 'DAMulti_SMAP',              # string used in output file names for sums & stats,
                                                       #   can be same as expid or different -- e.g., reflect info
                                                       #   about subset of included species or about cross-masking
        'domain'        : 'SMAP_EASEv2_M36_GLOBAL',  
        'da_t0'         : 3,                           # (fractional) UTC hour of first land analysis 
        'da_dt'         : 10800,                       # ObsFcstAna file interval in seconds
        'species_list'  : [5,6,7,8],                   # indices of species to be processed 
        'obsparam_time' : "20150401_0000"              # time stamp of obsparam file (YYYYMMDD_HHMM)
    }

    # Optional experiment(s) can be added for cross-masking or extracting obs from a different experiment.
    #
    # All optional experiments and the main experiment must have identical tile space (BCs resolution) and domain.
    #
    # If the default "species" number/order do not match, set "species_list" accordingly to force a match.
    # Output will be cross-masked between all specified experiments.    
    #
    # Forecasts and analyses are always from the main experiment.
    # Observations are from the experiment with 'use_obs' set to True (default is exp_main).  The most
    #   likely use case for reading obs from a supplemental experiment is when computing OmF etc diagnostics
    #   for an open loop experiment that only has unscaled obs, and _scaled_ obs must be read from a 
    #   coresponding DA experiment.

    exp_sup1 = {
        'expdir'        : '/discover/nobackup/projects/gmao/merra/iau/merra_land/SMAP_runs/SMAP_Nature_v11/',
        'expid'         : 'DAv8_M36',
        'use_obs'       : True,                        # if True, use obs data from this exp
        'species_list'  : [1,2,3,4],                   # indices of species to be processed,
                                                       #   must identify same species as selected in main exp
        'obsparam_time' : "20150401_0000"              # time stamp of obsparam file (YYYYMMDD_HHMM)
    }

    # Convert experiments input to a list; first entry must be exp_main: 

    #exp_list = [exp_main]               # no cross-masking
    exp_list = [exp_main, exp_sup1]      # cross-mask exp_main with exp_sup1
    
    # Top level directory for all output from this package:

    out_path = '/discover/nobackup/qliu/SMAP_test/'

    # Directory for monthly sum files:
    # - Can use the experiment directory or a different path.
    # - Automatically appends /Yyyyy/Mmm/ for individual months.
    
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

    domain = exp_list[0]['domain']
    
    for exp in exp_list:
        expdir        = exp['expdir']
        expid         = exp['expid']

        YYYY          = exp['obsparam_time'][0:4]
        MM            = exp['obsparam_time'][4:6]
        
        fop           = expdir+expid+'/output/'+domain+'/rc_out/Y'+YYYY+'/M'+MM+'/'+expid+'.ldas_obsparam.'+exp['obsparam_time']+'z.txt'
        obsparam_orig = read_obs_param(fop)

        # get the species list, default is all species 
        species_list = exp.get('species_list', [ int(obsparam_orig[i]['species']) for i in range(len(obsparam_orig)) ])
        
        # subset obsparam_orig based on species_list; keep order of obsparam_orig (independent of order of species_list)
        obsparam = []
        for i in range(len(obsparam_orig)):
            if int(obsparam_orig[i]['species']) in species_list:
                obsparam.append(obsparam_orig[i])
        
        ftc = expdir+expid+'/output/'+ domain+'/rc_out/'+ expid+'.ldas_tilecoord.bin'
        tc  = read_tilecoord(ftc)

        ftg = expdir+expid+'/output/'+ domain+'/rc_out/'+ expid+'.ldas_tilegrids.bin'
        tg_global, tg_domain  = read_tilegrids(ftg)

        # add tilecoord and obsparam into to exp        
        exp.update({'tilecoord':tc, 'obsparam':obsparam, 'tilegrid_global':tg_global,'tilegrid_domain': tg_domain})

    # verify that obs species match across experiments
    for exp_idx, exp in enumerate(exp_list) :
        obsparam = exp.get('obsparam')
        if exp_idx==0 :
            obsparam_main = obsparam
        else:
            if len(obsparam) != len(obsparam_main) :
                print("ERROR: 'obsparam' mismatch (length)!  Check 'select_species' input in user_config.py." )
                sys.exit()
            else:
                for a, b in zip(obsparam, obsparam_main) :
                    if a['descr'] != b['descr'] :
                        print("ERROR: 'obsparam' mismatch (descr)!  Check 'select_species' input in user_config.py." )
                        sys.exit()

    # wrap up config
    config ={
        'exp_list'   : exp_list,
        'start_time' : start_time,
        'end_time'   : end_time,
        'sum_path'   : sum_path,
        'out_path'   : out_path,
        }

    return config

# ====================== EOF =========================================================
