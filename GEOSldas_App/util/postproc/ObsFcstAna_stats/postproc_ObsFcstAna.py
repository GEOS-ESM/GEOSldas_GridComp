
# Tool to prepocess GEOSldas ObsFcstAna output into monthly sums, sums of squares, and sum of cross-products
#
# qliu, amfox, rreichle - May 2025

import numpy as np
import os
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from netCDF4 import  Dataset, date2num
import sys

sys.path.append('../../shared/python/')

from util                         import make_folder
from helper.compute_monthly_stats import compute_monthly_stats
from helper.write_nc4             import write_sums_nc4, write_stats_nc4, write_omf_stats_nc4, write_omf_grouped_stats_nc4

import warnings; warnings.filterwarnings("ignore")

class postproc_ObsFcstAna:
    
    def __init__(self, exp_list, start_time, end_time, obs_from=0):
        self.expdir_list   = [item['expdir'] for item in exp_list]
        self.expid_list    = [item['expid']  for item in exp_list]
        self.exptag_list   = [item['exptag'] for item in exp_list]
        self.domain        = exp_list[0]['domain']
        self.start_time    = start_time 
        self.end_time      = end_time 
        self.var_list      = ['obs_obs','obs_obsvar','obs_fcst','obs_fcstvar','obs_ana','obs_anavar']
        self.tilecoord     = exp_list[0]['tilecoord']
        self.obsparam_list = [item['obsparam'] for item in exp_list]
        self.obs_from      = obs_from

    # ---------------------------------------------------------------------------
        
    def save_monthly_sum(self, outpath='./'):
        expdir_list   = self.expdir_list
        expid_list    = self.expid_list
        exptag_list   = self.exptag_list
        domain        = self.domain
        var_list      = self.var_list
        tc            = self.tilecoord
        obsparam_list = self.obsparam_list
        start_time    = self.start_time
        end_time      = self.end_time

        start_month =  start_time.replace(day=1,hour=3)
        if end_time.day ==1:
            end_month =  end_time.replace(day=1,hour=3)
        else:
            tmp_time = end_time + relativedelta(months =1)
            end_month =tmp_time.replace(day=1,hour=3)
        
        current_month = start_month
        # month loop
        while current_month < end_month:
            
            # make monthly file output directory 
            mo_path = outpath + '/Y'+ current_month.strftime('%Y') + '/M' + current_month.strftime('%m') + '/'
            make_folder(mo_path)
            
            outid =  '_'.join([item for item in exptag_list])

            if self.obs_from  > 0:
                outid = outid + '_Obs_from_' + exptag_list[self.obs_from]
                
            if current_month == start_month and start_time.day > 1:
                sdate = start_time
                edate = current_month + relativedelta(months =1)
                fout = outid + '.ens_avg.ldas_ObsFcstAna.from_' + start_time.strftime('%Y%m%d') +'_stats.nc4'
            elif current_month.year == end_time.year and current_month.month == end_time.month:
                sdate  = current_month
                edate = end_time
                fout = outid+ '.ens_avg.ldas_ObsFcstAna.to_' + end_time.strftime('%Y%m%d') +'_stats.nc4'
            else:
                sdate = current_month
                edate = current_month + relativedelta(months =1)
                fout = outid + '.ens_avg.ldas_ObsFcstAna.' + current_month.strftime('%Y%m') +'_stats.nc4'

            fout = mo_path + fout
            
            # allow stats based on partial month
            month_range = [sdate, edate]
            # skip if output file already exists
            if  not os.path.isfile(fout):
                print('computing monthly sums ...')
                # compute monthly sum
                mN_data, mdata_sum, mdata2_sum, moxf_sum, moxa_sum, mfxa_sum = \
                        compute_monthly_stats(expdir_list,expid_list,domain,month_range,tc,obsparam_list,var_list,self.obs_from)

                # save monthly sum in nc4 file
                write_sums_nc4(fout, mN_data,mdata_sum, mdata2_sum, moxf_sum, moxa_sum, mfxa_sum, obsparam_list[0])
            else:
                print('file exist, skip '+fout)
                
            current_month = current_month + relativedelta(months=1)      

    # ---------------------------------------------------------------------------
            
    def calculate_stats_fromsums(self, mo_path='./', write_to_nc=True, filename='./stats.nc4'):
        
        start_time  = self.start_time
        end_time    = self.end_time
        exptag_list = self.exptag_list
        
        # Variable list for computing sum and sum of squared
        var_list = self.var_list 

        # Read tilecoord and obsparam for tile and obs species information
        n_tile = self.tilecoord['N_tile']
        n_spec = len(self.obsparam_list[0])

        # Initialize statistical metrics 
        data_sum  = {}
        data2_sum = {}
        N_data  = np.zeros((n_tile, n_spec))
        oxf_sum = np.zeros((n_tile, n_spec))
        oxa_sum = np.zeros((n_tile, n_spec))
        fxa_sum = np.zeros((n_tile, n_spec))

        for var in var_list:
            data_sum[var] = np.zeros((n_tile, n_spec))
            data2_sum[var] = np.zeros((n_tile, n_spec))

        # Time loop: processing data at monthly time step
        
        start_month =  start_time.replace(day=1,hour=3)
        if end_time.day ==1:
            end_month =  end_time.replace(day=1,hour=3)
        else:
            tmp_time = end_time + relativedelta(months =1)
            end_month =tmp_time.replace(day=1,hour=3)
        
        current_month = start_month
        # month loop
        while current_month < end_month:
            
            fpath = mo_path + '/Y'+ current_month.strftime('%Y') + '/M' + current_month.strftime('%m') + '/'

            outid =  '_'.join([item for item in exptag_list])

            if self.obs_from  > 0:
                outid = outid + '_Obs_from_' + exptag_list[self.obs_from]
                
            if current_month == start_month and start_time.day > 1:
                sdate = start_time
                edate = current_month + relativedelta(months =1)
                fout = outid + '.ens_avg.ldas_ObsFcstAna.from_' + start_time.strftime('%Y%m%d') +'_stats.nc4'
            elif current_month.year == end_time.year and current_month.month == end_time.month:
                sdate  = current_month
                edate = end_time
                fout = outid+ '.ens_avg.ldas_ObsFcstAna.to_' + end_time.strftime('%Y%m%d') +'_stats.nc4'
            else:
                sdate = current_month
                edate = current_month + relativedelta(months =1)
                fout = outid + '.ens_avg.ldas_ObsFcstAna.' + current_month.strftime('%Y%m') +'_stats.nc4'

            fout = fpath + fout
                            
            # Read monthly data if file exists, otherwise compute monthly statistics first   
            if os.path.isfile(fout):
                print('read sums from  monthly file: '+fout)
                mdata_sum = {}
                mdata2_sum = {}
                with Dataset(fout,'r') as nc:
                    mN_data  = nc.variables['N_data'][:]
                    moxf_sum = nc.variables['obsxfcst_sum'][:]
                    moxa_sum = nc.variables['obsxana_sum'][:]
                    mfxa_sum = nc.variables['fcstxana_sum'][:]
                    for var in var_list:
                        mdata_sum[var] = nc.variables[var+'_sum'][:]
                        mdata2_sum[var] = nc.variables[var+'2_sum'][:]
                       
                # Aggregate monthly data
                N_data  += mN_data
                oxf_sum += moxf_sum
                oxa_sum += moxa_sum
                fxa_sum += mfxa_sum
               
                for var in var_list:
                    data_sum[var] += mdata_sum[var] 
                    data2_sum[var] += mdata2_sum[var]  
            else:
                raise FileNotFoundError(f"File {fout} does not exist")
                 
            current_month =current_month + relativedelta(months=1)

        # Compute the basic statistics after finishing time loop based on the accumulated data.
        data_mean  = {}
        data2_mean = {}
        data_var   = {}

        # calculation  
        for var in var_list:
            data_sum[var][ N_data == 0] = np.nan
            data2_sum[var][N_data == 0] = np.nan
            
            data_mean[ var] = data_sum[var]  / N_data
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

        stats = {}
        for var in var_list:
            stats[var[4:]+'_mean']     = data_mean[var]
            stats[var[4:]+'_variance'] = data_var[ var]
        stats['oxf_mean'] = oxf_mean
        stats['oxa_mean'] = oxa_mean
        stats['fxa_mean'] = fxa_mean
        stats['N_data']   = N_data

        if write_to_nc:
            print('writing stats nc4 file: '+filename)
            write_stats_nc4(filename, stats)
            
        return stats

    # ---------------------------------------------------------------------------
    
    def calculate_monthly_omf(self, mo_path='./', write_to_nc=True, filename='./omf_stats.nc4'):
        """
        Compute monthly observation-minus-forecast (OmF) and observation-minus-analysis (OmA) 
        statistics from GEOSldas ObsFcstAna diagnostics. 
        
        For each month, the function reads precomputed sum and sum-of-squares fields, 
        calculates per-tile means, variances, covariances, and diagnostic metrics, 
        including normalized OmF statistics. Optionally writes results to NetCDF.

        Parameters:
        -----------
        mo_path : str
            Path to monthly summary files (organized as mo_path/Y<year>/M<month>/...).
        
        write_to_nc : bool
            If True, writes the aggregated statistics to a NetCDF file.
        
        filename : str
            Output NetCDF filename if write_to_nc is True.

        Returns:
        --------
        omf_stats_combined : dict
            Dictionary containing monthly time series (stacked by month) for each metric.
            Keys include:
                - 'N_data'            : Number of observations per tile
                - 'OmF_mean'          : Mean(Obs - Forecast)
                - 'OmF_stdv'          : Std(Obs - Forecast)
                - 'OmF_norm_mean'     : Normalized OmF mean
                - 'OmF_norm_stdv'     : Normalized OmF std
                - 'OmA_mean'          : Mean(Obs - Analysis)
                - 'OmA_stdv'          : Std(Obs - Analysis)
                - 'time'              : List of datetime objects for each month
        """

        start_time  = self.start_time
        end_time    = self.end_time
        exptag_list = self.exptag_list
        var_list    = self.var_list  # variables to sum/sum-of-squares

        Nmin = 5

        # Prepare storage for combined monthly diagnostics
        month_list = []
        combined = {
            'N_data': [], 'OmF_mean': [], 'OmF_stdv': [],
            'OmF_norm_mean': [], 'OmF_norm_stdv': [],
            'OmA_mean': [], 'OmA_stdv': []
        }

        # Determine month loop bounds
        start_month = start_time.replace(day=1, hour=3)
        if end_time.day == 1:
            end_month = end_time.replace(day=1, hour=3)
        else:
            end_month = (end_time + relativedelta(months=1)).replace(day=1, hour=3)
        current_month = start_month

        # Loop over each month
        while current_month < end_month:
            # Build file path and output name for sums
            fpath = os.path.join(mo_path,
                                'Y' + current_month.strftime('%Y'),
                                'M' + current_month.strftime('%m'))
            outid = '_'.join(exptag_list)
            if self.obs_from > 0:
                outid += '_Obs_from_' + exptag_list[self.obs_from]

            # Set monthly filename
            if current_month == start_month and start_time.day > 1:
                fout_name = f"{outid}.ens_avg.ldas_ObsFcstAna.from_{start_time.strftime('%Y%m%d')}_stats.nc4"
            elif (current_month.year == end_time.year and
                current_month.month == end_time.month):
                fout_name = f"{outid}.ens_avg.ldas_ObsFcstAna.to_{end_time.strftime('%Y%m%d')}_stats.nc4"
            else:
                fout_name = f"{outid}.ens_avg.ldas_ObsFcstAna.{current_month.strftime('%Y%m')}_stats.nc4"

            fout = os.path.join(fpath, fout_name)

            # Read this month's sums
            if not os.path.isfile(fout):
                raise FileNotFoundError(f"File {fout} does not exist")
            print('read sums from monthly file:', fout)

            # Load monthly accumulators
            with Dataset(fout, 'r') as nc:
                N_data    = nc.variables['N_data'][:]
                oxf_sum   = nc.variables['obsxfcst_sum'][:]
                oxa_sum   = nc.variables['obsxana_sum'][:]
                fxa_sum   = nc.variables['fcstxana_sum'][:]
                data_sum  = {var: nc.variables[var + '_sum'][:] for var in var_list}
                data2_sum = {var: nc.variables[var + '2_sum'][:] for var in var_list}

            # Compute the basic statistics after finishing time loop based on the accumulated data.
            data_mean ={}
            data2_mean = {}
            data_var = {}

            # calculation  
            for var in var_list:
                data_sum[var][ N_data == 0] = np.nan
                data2_sum[var][N_data == 0] = np.nan

                data_mean[var]  = data_sum[var]  / N_data
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

            stats = {}
            for var in var_list:
                stats[var[4:]+'_mean']= data_mean[var]
                stats[var[4:]+'_variance'] = data_var[var]
            stats['oxf_mean'] = oxf_mean
            stats['oxa_mean'] = oxa_mean
            stats['fxa_mean'] = fxa_mean
            stats['N_data']   = N_data

            # Then computer metrics of O-F, O-A, etc. based on above computed
            N_data = stats['N_data']

            # mean(x-y) = E[x] - E[y]   
            OmF_mean = stats['obs_mean'] - stats['fcst_mean']
            OmA_mean = stats['obs_mean'] - stats['ana_mean']
            # var(x-y) = var(x) + var(y) - 2cov(x,y)
            # cov(x,y) = E[xy] - E[x]E[y]
            OmF_stdv  = np.sqrt(np.maximum(0, stats['obs_variance'] + stats['fcst_variance'] - \
                                2 * (stats['oxf_mean'] - stats['obs_mean']*stats['fcst_mean'])))

            OmA_stdv  = np.sqrt(np.maximum(0, stats['obs_variance'] + stats['ana_variance'] - \
                                2 * (stats['oxa_mean'] - stats['obs_mean']*stats['ana_mean'])))

            # "fcstvar" is assumed constant here for convenience. Modify if necessary
            OmF_norm_mean = OmF_mean / np.sqrt(stats['obsvar_mean'] + stats['fcstvar_mean']) 
            OmF_norm_stdv = np.sqrt(OmF_stdv**2 / (stats['obsvar_mean'] + stats['fcstvar_mean']) )

            # Mask out data points with insufficent observations using the Nmin threshold
            # Do NOT apply to N_data
            OmF_mean[     N_data < Nmin] = np.nan
            OmF_stdv[     N_data < Nmin] = np.nan
            OmF_norm_mean[N_data < Nmin] = np.nan
            OmF_norm_stdv[N_data < Nmin] = np.nan
            OmA_mean[     N_data < Nmin] = np.nan
            OmA_stdv[     N_data < Nmin] = np.nan

            omf_stats = {
                'N_data':        N_data,
                'OmF_mean':      OmF_mean,
                'OmF_stdv':      OmF_stdv,
                'OmF_norm_mean': OmF_norm_mean,
                'OmF_norm_stdv': OmF_norm_stdv,
                'OmA_mean':      OmA_mean,
                'OmA_stdv':      OmA_stdv
            }

            # Store into combined container
            for key in combined:
                combined[key].append(omf_stats[key].copy())
            month_list.append(current_month)

            # Increment to the next month
            current_month = current_month + relativedelta(months=1)

        # Stack combined stats into arrays
        omf_stats_combined = {k: np.stack(v, axis=0) for k, v in combined.items()}
        omf_stats_combined['time'] = np.array(month_list)

        if write_to_nc:
            print('writing omf stats nc4 file: '+filename)
            write_omf_stats_nc4(filename, omf_stats_combined)

        return omf_stats_combined

    # ---------------------------------------------------------------------------
    
    def calculate_monthly_omf_by_sensor(self, sensor_groups, mo_path='./', write_to_nc=True, filename='./omf_group_stats.nc4'):
        
        """
        Compute per-tile monthly OmF and OmA statistics grouped by sensor.

        Parameters:
            sensor_groups : dict
                Mapping from group name (e.g., 'SMAP') to list of species indices.
            mo_path : str
                Path to monthly input files.
            write_to_nc : bool
                Whether to write results to NetCDF.
            filename : str
                Path to output NetCDF file.

        Returns:
            group_stats : dict
                Dictionary of per-group, per-tile monthly statistics, with shape (time, tile).
        """        
       
        # Setup
        start_time  = self.start_time
        end_time    = self.end_time
        exptag_list = self.exptag_list
        var_list    = self.var_list

        Nmin = 10

        # Time bounds
        start_month = start_time.replace(day=1, hour=3)
        end_month = (end_time + relativedelta(months=1)).replace(day=1, hour=3) if end_time.day != 1 else end_time.replace(day=1, hour=3)

        # Initialize storage
        month_list = []
        group_combined = {
            g: {key: [] for key in ['N_data', 'OmF_mean', 'OmF_stdv', 'OmF_norm_mean', 'OmF_norm_stdv', 'OmA_mean', 'OmA_stdv']}
            for g in sensor_groups
        }

        current_month = start_month
        while current_month < end_month:
            # File path
            fpath = os.path.join(mo_path, 'Y' + current_month.strftime('%Y'), 'M' + current_month.strftime('%m'))
            outid = '_'.join(exptag_list)
            if self.obs_from > 0:
                outid += '_Obs_from_' + exptag_list[self.obs_from]

            if current_month == start_month and start_time.day > 1:
                fname = f"{outid}.ens_avg.ldas_ObsFcstAna.from_{start_time.strftime('%Y%m%d')}_stats.nc4"
            elif current_month.year == end_time.year and current_month.month == end_time.month:
                fname = f"{outid}.ens_avg.ldas_ObsFcstAna.to_{end_time.strftime('%Y%m%d')}_stats.nc4"
            else:
                fname = f"{outid}.ens_avg.ldas_ObsFcstAna.{current_month.strftime('%Y%m')}_stats.nc4"

            fout = os.path.join(fpath, fname)
            if not os.path.isfile(fout):
                raise FileNotFoundError(f"Missing file: {fout}")
            print('Reading:', fout)

            # Load monthly accumulators
            with Dataset(fout, 'r') as nc:
                N_data    = nc.variables['N_data'][:]
                oxf_sum   = nc.variables['obsxfcst_sum'][:]
                oxa_sum   = nc.variables['obsxana_sum'][:]
                data_sum  = {var: nc.variables[var + '_sum'][:] for var in var_list}
                data2_sum = {var: nc.variables[var + '2_sum'][:] for var in var_list}

            # --- Group-level analysis only ---
            for group, indices in sensor_groups.items():
                N_group = np.sum(N_data[:, indices], axis=1).astype(float)
                N_group[N_group == 0] = np.nan

                # Mean and variance
                mean  = {var[4:]: np.sum(data_sum[var][:, indices], axis=1) / N_group for var in var_list}
                mean2 = {var[4:]: np.sum(data2_sum[var][:, indices], axis=1) / N_group for var in var_list}
                var_  = {var[4:]: mean2[var[4:]] - mean[var[4:]]**2 for var in var_list}

                # Covariance terms
                oxf_mean = np.sum(oxf_sum[:, indices], axis=1) / N_group
                oxa_mean = np.sum(oxa_sum[:, indices], axis=1) / N_group

                OmF_mean = mean['obs'] - mean['fcst']
                OmF_stdv = np.sqrt(np.maximum(0, var_['obs'] + var_['fcst'] - 2 * (oxf_mean - mean['obs'] * mean['fcst'])))

                OmA_mean = mean['obs'] - mean['ana']
                OmA_stdv = np.sqrt(np.maximum(0, var_['obs'] + var_['ana'] - 2 * (oxa_mean - mean['obs'] * mean['ana'])))

                # "fcstvar" is assumed constant here for convenience. Modify if necessary
                OmF_norm_mean = OmF_mean / np.sqrt(mean['obsvar'] + mean['fcstvar'])
                OmF_norm_stdv = np.sqrt(OmF_stdv**2 / (mean['obsvar'] + mean['fcstvar']) )                          

                # Mask out data points with insufficent observations using the Nmin threshold
                # Do NOT apply to N_data
                OmF_mean[     N_group < Nmin] = np.nan
                OmF_stdv[     N_group < Nmin] = np.nan
                OmA_mean[     N_group < Nmin] = np.nan
                OmA_stdv[     N_group < Nmin] = np.nan
                OmF_norm_mean[N_group < Nmin] = np.nan
                OmF_norm_stdv[N_group < Nmin] = np.nan

                # Aggregate monthly data
                group_combined[group]['N_data'       ].append(N_group)
                group_combined[group]['OmF_mean'     ].append(OmF_mean)
                group_combined[group]['OmF_stdv'     ].append(OmF_stdv)
                group_combined[group]['OmF_norm_mean'].append(OmF_norm_mean)
                group_combined[group]['OmF_norm_stdv'].append(OmF_norm_stdv)
                group_combined[group]['OmA_mean'     ].append(OmA_mean)
                group_combined[group]['OmA_stdv'     ].append(OmA_stdv)

            month_list.append(current_month)
            current_month += relativedelta(months=1)

        # Final stack
        group_stats = {
            group: {k: np.stack(v, axis=0) for k, v in stats.items()}
            for group, stats in group_combined.items()
        }
        group_stats['time'] = np.array(month_list)


        if write_to_nc:
            print("Writing netCDF:", filename)
            write_omf_grouped_stats_nc4(filename, group_stats)

        return group_stats

# ============== EOF ====================================================================
