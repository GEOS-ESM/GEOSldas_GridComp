
# Tool to preprocess GEOSldas ObsFcstAna output into monthly sums, sums of squares, and sums of cross-products
#
# qliu, amfox, rreichle - May 2025

import numpy as np
import os
import yaml

import sys;       sys.path.append('../../shared/python/')

import warnings;  warnings.filterwarnings("ignore")

from datetime                    import datetime, timedelta
from dateutil.relativedelta      import relativedelta
from netCDF4                     import Dataset, date2num
from read_GEOSldas               import read_ObsFcstAna, read_tilecoord, read_obs_param

from helper.write_nc4            import write_sums_nc4, write_stats_nc4

class postproc_ObsFcstAna:
    
    def __init__(self, exp_list, start_time, end_time, sum_path='./'):
        self.exp_list        = exp_list
        self.expdir_list     = [item['expdir'] for item in exp_list]
        self.expid_list      = [item['expid']  for item in exp_list]
        self.exptag          = exp_list[0]['exptag']
        self.domain          = exp_list[0]['domain']
        self.start_time      = start_time 
        self.end_time        = end_time
        self.da_t0           = exp_list[0]['da_t0']
        self.da_dt           = exp_list[0]['da_dt']
        self.var_list        = ['obs_obs','obs_obsvar','obs_fcst','obs_fcstvar','obs_ana','obs_anavar']
        self.tilecoord       = exp_list[0]['tilecoord']
        self.tilegrid_global = exp_list[0]['tilegrid_global']
        self.tilegrid_domain = exp_list[0]['tilegrid_domain']
        self.obsparam_list   = [item['obsparam'] for item in exp_list]
        self.sum_path        = sum_path

        # Determine experiment that supplies obs data
        self.obs_from = -1
        for exp_idx, exp in enumerate(exp_list):
            if exp.get('use_obs', None):              # found use_obs=True 
                if self.obs_from >= 0:          
                    print("ERROR: use_obs=True in multiple experiments. Edit user_config.py to remove conflict." )
                    sys.exit()
                else:
                    self.obs_from = exp_idx
        if self.obs_from < 0: self.obs_from = 0       # by default, obs data are from exp_list[0]
        print(f"obs data are from {exp_list[self.obs_from]['expid']}")   

        # Verify the configuration every time when current class is initialized
        # to avoid saving sums with different configs in the same directory

        # Same configurations should have identical values for these fields
        config_verify = ['expdir','expid','exptag','domain','use_obs','species_list']

        # Construct config for each experiment
        config_list = []
        for exp in exp_list:
            config_list.append({var:exp[var] for var in config_verify if var in exp})

        # File of configuration for verification
        f_config = self.sum_path + '/' + self.exptag + '_config.yaml'

        # Save a new file or compare current configuration with previously saved 
        if not os.path.exists(f_config):
            with open(f_config, 'w') as f:
                yaml.dump(config_list, f, default_flow_style=False)
            print(f'Configuration saved to {f_config}')
        else:
            with open(f_config,'r') as f:
                saved_exp_list = yaml.safe_load(f)
            print(f'Found saved configuration')

            if config_list != saved_exp_list:
                print('user configuration is different from previously saved '+f_config)
                sys.exit()  
    # ----------------------------------------------------------------------------------------------------------
    #
    # Function to compute monthly sums from x-hourly ObsFcstAna data for a single month.

    def compute_monthly_sums(self, date_time):
        
        expdir_list   = self.expdir_list
        expid_list    = self.expid_list
        tc            = self.tilecoord
        obsparam_list = self.obsparam_list
        var_list      = self.var_list
        obs_from      = self.obs_from
        da_dt         = self.da_dt
        
        n_tile = tc['N_tile']
        n_spec = len(obsparam_list[0])

        date_time = date_time.replace(hour=self.da_t0)
        end_time  = date_time + relativedelta(months=1) 
        
        data_sum  = {}
        data2_sum = {}

        N_data  = np.zeros((n_tile, n_spec))
        oxf_sum = np.zeros((n_tile, n_spec))
        oxa_sum = np.zeros((n_tile, n_spec))
        fxa_sum = np.zeros((n_tile, n_spec))

        for var in var_list:
            data_sum[ var] = np.zeros((n_tile, n_spec))
            data2_sum[var] = np.zeros((n_tile, n_spec))

        while date_time < end_time:
            
            # read the list of experiments at each time step (OFA="ObsFcstAna")
            OFA_list = []
            for i in range(len(expdir_list)):
                fname = expdir_list[i]+expid_list[i]+'/output/'+self.domain+'/ana/ens_avg/Y'+ \
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
                    data_tile[var] = np.zeros((n_tile, n_spec)) + np.nan

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

            # cross-mask over all experiments 
            is_cross_valid = ~np.isnan(data_all[0]['obs_obs'])
            for data in data_all[1:]:
                mask = ~np.isnan(data['obs_obs'])
                is_cross_valid = np.logical_and(is_cross_valid,mask)

            # reconstruct the output variable dictionary based on input options;
            # obs_obs and obs_obsvar are from exp_list[obs_from], the rest are from exp_list[0]
            data_tile = {}    
            for var in var_list:
                if 'obs_obs' in var:
                    data_tile[var] = data_all[obs_from][var]
                else:
                    data_tile[var] = data_all[0][var]

            is_valid = is_cross_valid
            
            N_data[ is_valid] += 1
            oxf_sum[is_valid] += data_tile['obs_obs' ][is_valid] * data_tile['obs_fcst'][is_valid]
            oxa_sum[is_valid] += data_tile['obs_obs' ][is_valid] * data_tile['obs_ana' ][is_valid]
            fxa_sum[is_valid] += data_tile['obs_fcst'][is_valid] * data_tile['obs_ana' ][is_valid]

            for var in var_list:
                data_sum[ var][is_valid] += data_tile[var][is_valid]
                data2_sum[var][is_valid] += data_tile[var][is_valid] **2
            
            date_time = date_time + timedelta(seconds=da_dt)

        return N_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum

    # ----------------------------------------------------------------------------------------------------------
    #
    # Function to compute monthly sums and save results in nc4 files for all months in [start_time, end_time].
    #
    # Skips computation/saving of monthly sums output if file already exists.

    def save_monthly_sums(self):
        expdir_list    = self.expdir_list
        expid_list     = self.expid_list
        var_list       = self.var_list
        tc             = self.tilecoord
        obsparam_list  = self.obsparam_list

        date_time   = self.start_time

        while date_time < self.end_time:           # loop through months
            
            # make monthly file output directory 
            mo_path = self.sum_path + '/Y'+ date_time.strftime('%Y') + '/M' + date_time.strftime('%m') + '/'
            os.makedirs(mo_path, exist_ok=True)
    
            fout = self.exptag + '.ens_avg.ldas_ObsFcstAna_sums.' + date_time.strftime('%Y%m') +'.nc4'

            fout = mo_path + fout
            
            # skip if output file already exists
            if  not os.path.isfile(fout):
                print('computing monthly sums...')
                # compute monthly sums
                mN_data, mdata_sum, mdata2_sum, moxf_sum, moxa_sum, mfxa_sum = \
                    self.compute_monthly_sums(date_time)

                # save monthly sums in nc4 file
                write_sums_nc4(fout, mN_data,mdata_sum, mdata2_sum, moxf_sum, moxa_sum, mfxa_sum, obsparam_list[0])
            else:
                print('file exists, skipping '+fout)
                
            date_time = date_time + relativedelta(months=1)      

    # ----------------------------------------------------------------------------------------------------------
    #
    # Function to compute long-term temporal statistics of individual species based on monthly sums.
    #
    # Assumes that monthly sums files have been saved [see save_monthly_sums()].

    def calc_temporal_stats_from_sums(self, write_to_nc=True, fout_stats='./stats.nc4'):

        # Variable list for computing sum and sum-of-squares
        var_list = self.var_list 

        # Read tilecoord and obsparam for tile and obs species information
        n_tile = self.tilecoord['N_tile']
        n_spec = len(self.obsparam_list[0])

        # ---------------------------------------------------------------
        #
        # compute accumulated sums for period (start_time, end_time)
        
        # Initialize statistical metrics 
        data_sum  = {}
        data2_sum = {}
        N_data    = np.zeros((n_tile, n_spec))
        oxf_sum   = np.zeros((n_tile, n_spec))
        oxa_sum   = np.zeros((n_tile, n_spec))
        fxa_sum   = np.zeros((n_tile, n_spec))

        for var in var_list:
            data_sum[ var] = np.zeros((n_tile, n_spec))
            data2_sum[var] = np.zeros((n_tile, n_spec))

        # Time loop: processing data at monthly time step
        
        date_time = self.start_time

        while date_time < self.end_time:      # loop through months
            
            fpath = self.sum_path + '/Y'+ date_time.strftime('%Y') + '/M' + date_time.strftime('%m') + '/'
                
            fout  = self.exptag + '.ens_avg.ldas_ObsFcstAna_sums.' + date_time.strftime('%Y%m') +'.nc4'

            fout  = fpath + fout
                            
            # Read and accumulate monthly sums data 
            if os.path.isfile(fout):
                print('read sums from monthly file: '+fout)
                mdata_sum  = {}
                mdata2_sum = {}
                with Dataset(fout,'r') as nc:
                    mN_data             = nc.variables['N_data'      ][:]
                    moxf_sum            = nc.variables['obsxfcst_sum'][:]
                    moxa_sum            = nc.variables['obsxana_sum' ][:]
                    mfxa_sum            = nc.variables['fcstxana_sum'][:]
                    for var in var_list:
                        mdata_sum[ var] = nc.variables[var+'_sum'    ][:]
                        mdata2_sum[var] = nc.variables[var+'2_sum'   ][:]
                       
                # Aggregate monthly data
                N_data  += mN_data
                oxf_sum += moxf_sum
                oxa_sum += moxa_sum
                fxa_sum += mfxa_sum
               
                for var in var_list:
                    data_sum[ var] += mdata_sum[ var] 
                    data2_sum[var] += mdata2_sum[var]  
            else:
                raise FileNotFoundError(f"File {fout} does not exist, run save_monthly_sums() first")
                 
            date_time = date_time + relativedelta(months=1)

        # --------------------------------------------------------------------
        #
        # Compute stats (DA diagnostics) from accumulated sums data.
        
        data_mean  = {}
        data2_mean = {}
        data_var   = {}

        # calculation  
        for var in var_list:
            data_sum[var][ N_data == 0] = np.nan
            data2_sum[var][N_data == 0] = np.nan
            
            data_mean[ var] = data_sum[ var] / N_data
            data2_mean[var] = data2_sum[var] / N_data
            # var(x) = E[x2] - (E[x])^2
            data_var[var]   = data2_mean[var] - data_mean[var]**2
            
        oxf_sum[N_data == 0] = np.nan
        oxa_sum[N_data == 0] = np.nan
        fxa_sum[N_data == 0] = np.nan
        # E[xy]
        oxf_mean = oxf_sum / N_data
        oxa_mean = oxa_sum / N_data
        fxa_mean = fxa_sum / N_data

        # Then compute metrics of O-F, O-A, etc. based on above computed
        
        O_mean = data_mean['obs_obs']
        F_mean = data_mean['obs_fcst']
        A_mean = data_mean['obs_ana']

        O_stdv = np.sqrt(data_var['obs_obs'])
        F_stdv = np.sqrt(data_var['obs_fcst'])
        A_stdv = np.sqrt(data_var['obs_ana'])
        
        # mean(x-y) = E[x] - E[y]   
        OmF_mean = O_mean - F_mean
        OmA_mean = O_mean - A_mean
        
        # var(x-y) = var(x) + var(y) - 2cov(x,y)
        # cov(x,y) = E[xy] - E[x]E[y]
        OmF_stdv  = np.sqrt(O_stdv**2 + F_stdv**2 - 2 * (oxf_mean - O_mean*F_mean))
        OmA_stdv  = np.sqrt(O_stdv**2 + A_stdv**2 - 2 * (oxa_mean - O_mean*A_mean))

        # *****************************************************************************************
        # The time series mean and std-dev of the *normalized* OmF computed here are APPROXIMATED!
        # *****************************************************************************************
        # Here, we first compute the stats of the OmF time series and then normalize using 
        # the time-avg "obsvar" and "fcstvar" values.
        # Since "fcstvar" changes with time, the OmF values should be normalized at each time 
        # step (as in the older matlab scripts), and then the time series stats can be computed. 
        # To compute the exact stats with this python package, the sum and sum-of-squares of 
        # the normalized OmF values would need to be added into the sums files. 
        #
        OmF_norm_mean = OmF_mean / np.sqrt(data_mean['obs_obsvar'] + data_mean['obs_fcstvar'])      # APPROXIMATED stat!
        OmF_norm_stdv = np.sqrt(OmF_stdv**2 / (data_mean['obs_obsvar'] + data_mean['obs_fcstvar'])) # APPROXIMATED stat!
          
        # Mask out data points without any obs (NOTE: apply Nmin threshold when plotting or computing map avg)
        # Do NOT apply to N_data
        O_mean[       N_data < 1] = np.nan
        O_stdv[       N_data < 1] = np.nan
        F_mean[       N_data < 1] = np.nan
        F_stdv[       N_data < 1] = np.nan
        A_mean[       N_data < 1] = np.nan
        A_stdv[       N_data < 1] = np.nan
        
        OmF_mean[     N_data < 1] = np.nan
        OmF_stdv[     N_data < 1] = np.nan
        OmF_norm_mean[N_data < 1] = np.nan
        OmF_norm_stdv[N_data < 1] = np.nan
        OmA_mean[     N_data < 1] = np.nan
        OmA_stdv[     N_data < 1] = np.nan

        stats = {
            'O_mean'       : O_mean,        'O_stdv'       : O_stdv,
            'F_mean'       : F_mean,        'F_stdv'       : F_stdv,
            'A_mean'       : A_mean,        'A_stdv'       : A_stdv,
            'OmF_mean'     : OmF_mean,      'OmF_stdv'     : OmF_stdv,
            'OmA_mean'     : OmA_mean,      'OmA_stdv'     : OmA_stdv,
            'OmF_norm_mean': OmF_norm_mean, 'OmF_norm_stdv': OmF_norm_stdv,
            'N_data'       : N_data,
            }

        if write_to_nc:
            print('writing stats nc4 file: '+fout_stats)
            write_stats_nc4(fout_stats, stats)
            
        return stats

    # ----------------------------------------------------------------------------------------------------------
    #
    # Function to compute the O-F/O-A *spatial* statistics for a *single* month based on
    #   previously saved monthly sums.
    # Individual temporal and grid cell DA diagnostic values within a month are aggregated first;
    #   the monthly Ndata/mean/stdv are derived from the aggregated sample.  Consequently, in the
    #   final DA diagnostics, each obs value gets equal weight. Note that this differs from computing
    #   the straight spatial avg across a map of a given DA diagnostic.  The latter approach gives 
    #   the same weight to each location, regardless of how many obs are available at the location.

    def calc_spatial_stats_from_sums(self, date_time):

        var_list = ['obs_obs', 'obs_fcst','obs_ana']
        
        mo_path = self.sum_path + '/Y'+ date_time.strftime('%Y') + '/M' + date_time.strftime('%m') + '/'            
        fnc4_sums = mo_path + self.exptag + '.ens_avg.ldas_ObsFcstAna_sums.' + date_time.strftime('%Y%m') +'.nc4'
        
        mdata_sum = {}
        mdata2_sum = {}

        try:
            with Dataset(fnc4_sums,'r') as nc:
                mN_data                = nc.variables['N_data'      ][:]
                moxf_sum               = nc.variables['obsxfcst_sum'][:]
                moxa_sum               = nc.variables['obsxana_sum' ][:]
                moxf_sum[mN_data == 0] = np.nan
                moxa_sum[mN_data == 0] = np.nan
                for var in var_list:
                    mdata_sum[ var]               = nc.variables[var+'_sum' ][:]
                    mdata2_sum[var]               = nc.variables[var+'2_sum'][:]
                    mdata_sum[ var][mN_data == 0] = np.nan
                    mdata2_sum[var][mN_data == 0] = np.nan
        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.")
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)

        # Make sure only aggregate tiles with valid values for all variables
        for var in var_list:
            mN_data                           = mN_data.astype(float)
            mN_data[np.isnan(mdata_sum[var])] = np.nan
            mN_data[mN_data == 0]             = np.nan

        # cross mask before aggregating tile values   
        for var in var_list:
            mdata_sum[ var][np.isnan(mN_data)] = np.nan
            mdata2_sum[var][np.isnan(mN_data)] = np.nan
            moxf_sum[np.isnan(mN_data)]        = np.nan
            moxa_sum[np.isnan(mN_data)]        = np.nan

        # Aggregate data of all tiles
        N_data     = np.nansum(mN_data, axis=0)
        OxF_mean   = np.nansum(moxf_sum,axis=0)/N_data
        OxA_mean   = np.nansum(moxa_sum,axis=0)/N_data
        data_mean  = {}
        data2_mean = {}
        data_var   = {}
        for var in var_list:
            data_mean[ var] = np.nansum(mdata_sum[var ],axis=0)/N_data
            data2_mean[var] = np.nansum(mdata2_sum[var],axis=0)/N_data
            # var(x) = E[x2] - (E[x])^2
            data_var[var] = data2_mean[var] - data_mean[var]**2

        # Compute metrics of O-F, O-A, etc. based on above stats
        O_mean = data_mean['obs_obs']
        F_mean = data_mean['obs_fcst']
        A_mean = data_mean['obs_ana']

        O_var  = data_var[ 'obs_obs']
        F_var  = data_var[ 'obs_fcst']
        A_var  = data_var[ 'obs_ana']

        # mean(x-y) = E[x] - E[y]   
        OmF_mean = O_mean - F_mean
        OmA_mean = O_mean - A_mean

        # var(x-y) = var(x) + var(y) - 2cov(x,y)
        # cov(x,y) = E[xy] - E[x]E[y]
        OmF_stdv  = np.sqrt(O_var + F_var - 2 * (OxF_mean - O_mean*F_mean))
        OmA_stdv  = np.sqrt(O_var + A_var - 2 * (OxA_mean - O_mean*A_mean))

        stats = {
            'O_mean'  : O_mean,   'O_stdv'  : np.sqrt(O_var),
            'F_mean'  : F_mean,   'F_stdv'  : np.sqrt(F_var),
            'A_mean'  : A_mean,   'A_stdv'  : np.sqrt(A_var),
            'OmF_mean': OmF_mean, 'OmF_stdv': OmF_stdv,
            'OmA_mean': OmA_mean, 'OmA_stdv': OmA_stdv,
            'N_data': N_data,
            }

        return stats
    
# ============== EOF ====================================================================
