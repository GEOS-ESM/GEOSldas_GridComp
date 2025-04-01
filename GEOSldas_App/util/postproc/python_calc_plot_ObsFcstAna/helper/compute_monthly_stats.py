import numpy as np
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from helper.read_GEOSldas import read_ObsFcstAna, read_tilecoord, read_obs_param

def compute_monthly_stats(expdir,expid,domain,this_month,tc,obs_param,var_list):

    n_tile = tc['N_tile']
    n_spec = len(obs_param)

    start_time = this_month.replace(day=1,hour=3) 
    end_time = start_time + relativedelta(months=1)

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

        fname =  fname = expdir+expid+'/output/'+domain+'/ana/ens_avg/Y'+ \
                          date_time.strftime('%Y') + '/M' + \
                          date_time.strftime('%m') + '/'  + \
                          expid+'.ens_avg.ldas_ObsFcstAna.' + \
                          date_time.strftime('%Y%m%d_%H%M') +'z.bin'

        OFA = read_ObsFcstAna(fname)

        if len(OFA['obs_tilenum'] > 0):
            # Initialize full size variable to keep values
            data_tile={}
            for var in var_list:
                data_tile[var] = np.zeros((n_tile, n_spec)) +np.nan

            for ispec in np.arange(n_spec):
                # check species overall "assim" flag for masking
                this_species = int(obs_param[ispec]['species'])
                masked_data = {}
                if obs_param[ispec]['assim'] == 'T':
                    masked_tilenum = OFA['obs_tilenum'][np.logical_and(OFA['obs_species'] == this_species, OFA['obs_assim']==1)]
                    for var in var_list:
                        masked_data[var] = OFA[var][np.logical_and(OFA['obs_species'] == this_species, OFA['obs_assim']==1)]
                else:
                    masked_tilenum = OFA['obs_tilenum'][OFA['obs_species'] == this_species]
                    for var in var_list:
                        masked_data[var] = OFA[var][OFA['obs_species'] == this_species]          
     
                tile_idx = np.where(np.isin(tc['tile_id'], masked_tilenum))[0]

                for var in var_list:
                    data_tile[var][tile_idx, ispec] = masked_data[var]         

            is_valid = ~np.isnan(data_tile['obs_obs'])
            N_data[is_valid] += 1
            oxf_sum[is_valid] += data_tile['obs_obs'][is_valid] * data_tile['obs_fcst'][is_valid]
            oxa_sum[is_valid] += data_tile['obs_obs'][is_valid] * data_tile['obs_ana'][is_valid]
            fxa_sum[is_valid] += data_tile['obs_fcst'][is_valid] * data_tile['obs_ana'][is_valid]
            for var in var_list:
                data_sum[var][is_valid] += data_tile[var][is_valid]
                data2_sum[var][is_valid] += data_tile[var][is_valid] **2
        
        date_time = date_time + timedelta(seconds=10800)

    return N_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum

if __name__ == '__main__':
    date_time = datetime(2015,5,1)
    expdir = '/gpfsm/dnb05/projects/p51/SMAP_Nature/'
    expid = 'SPL4SM_Vv8010'
    domain = 'SMAP_EASEv2_M09_GLOBAL'
    var_list = ['obs_obs', 'obs_obsvar', 'obs_fcst', 'obs_fcstvar', 'obs_ana', 'obs_anavar']
    ftc = expdir+expid+'/output/'+domain+'/rc_out/'+expid+'.ldas_tilecoord.bin'
    tc = read_tilecoord(ftc)

    fop = expdir+expid+'/output/'+domain+'/rc_out/Y2015/M04/'+expid+'.ldas_obsparam.20150401_0000z.txt'
    obs_param = read_obs_param(fop)

    N_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum = \
           compute_monthly_stats(expdir,expid,domain,date_time,tc,obs_param,var_list)
