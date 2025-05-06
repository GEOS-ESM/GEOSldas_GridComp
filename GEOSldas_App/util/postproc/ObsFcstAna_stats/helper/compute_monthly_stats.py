import numpy as np
import sys
import os
sys.path.append('../../../shared/python/')

from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from read_GEOSldas import read_ObsFcstAna, read_tilecoord, read_obs_param

def compute_monthly_stats(expdir_list,expid_list,domain,this_month,tc,obsparam_list,var_list,obs_from):

    n_tile = tc['N_tile']
    n_spec = len(obsparam_list[0])

    #this_month = [start_date_ofthemonth, end_date_ofthemonth]
    start_time = this_month[0]
    end_time = this_month[1]
    
    data_sum = {}
    data2_sum = {}

    N_data = np.zeros((n_tile, n_spec))
    oxf_sum = np.zeros((n_tile, n_spec))
    oxa_sum = np.zeros((n_tile, n_spec))
    fxa_sum = np.zeros((n_tile, n_spec))

    for var in var_list:
        data_sum[var] = np.zeros((n_tile, n_spec))
        data2_sum[var] = np.zeros((n_tile, n_spec))

    date_time = start_time
    while date_time < end_time:
        
        # read the list of experiments at each time step
        OFA_list =[]
        for i in range(len(expdir_list)):
            fname = expdir_list[i]+expid_list[i]+'/output/'+domain+'/ana/ens_avg/Y'+ \
                              date_time.strftime('%Y') + '/M' + \
                              date_time.strftime('%m') + '/'  + \
                              expid_list[i]+'.ens_avg.ldas_ObsFcstAna.' + \
                              date_time.strftime('%Y%m%d_%H%M') +'z.bin'

            if os.path.isfile(fname):
                print('read '+fname)
                OFA_list.append(read_ObsFcstAna(fname))

        data_all=[]
        for OFA, obs_param in zip(OFA_list,obsparam_list):

            # Initialize full size variable for an experiment
            data_tile={}
            for var in var_list:
                data_tile[var] = np.zeros((n_tile, n_spec)) +np.nan

            if len(OFA['obs_tilenum']) > 0:
                for ispec in np.arange(n_spec):
                    # check species overall "assim" flag for masking
                    this_species = int(obs_param[ispec]['species'])
                    masked_data = {}

                    # get mask 
                    if obs_param[ispec]['assim'] == 'T':
                        is_valid = np.logical_and(OFA['obs_species'] == this_species, OFA['obs_assim']==1)
                    else:
                        is_valid = OFA['obs_species'] == this_species

                    tile_idx = OFA['obs_tilenum'][is_valid]-1
                    for var in var_list:
                        masked_data[var] = OFA[var][is_valid]
          
                    for var in var_list:
                        data_tile[var][tile_idx, ispec] = masked_data[var]         

            data_all.append(data_tile)

        # cross mask over all experiments 
        is_cross_valid = ~np.isnan(data_all[0]['obs_obs'])
        for data in data_all[1:]:
            mask = ~np.isnan(data['obs_obs'])
            is_cross_valid = np.logical_and(is_cross_valid,mask)

        # reconstruct the output variable dictionary based on input options
        # obs_obs and obs_obsvar are from exp_list[obs_from], the rest are from exp_list[0]
        data_tile = {}    
        for var in var_list:
            if 'obs_obs' in var:
                data_tile[var] = data_all[obs_from][var]
            else:
                data_tile[var] = data_all[0][var]

        is_valid = is_cross_valid       
        N_data[is_valid] += 1
        oxf_sum[is_valid] += data_tile['obs_obs'][is_valid] * data_tile['obs_fcst'][is_valid]
        oxa_sum[is_valid] += data_tile['obs_obs'][is_valid] * data_tile['obs_ana'][is_valid]
        fxa_sum[is_valid] += data_tile['obs_fcst'][is_valid] * data_tile['obs_ana'][is_valid]
        for var in var_list:
            data_sum[var][is_valid] += data_tile[var][is_valid]
            data2_sum[var][is_valid] += data_tile[var][is_valid] **2
        
        date_time = date_time + timedelta(seconds=10800)

    return N_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum

