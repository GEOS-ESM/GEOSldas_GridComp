import numpy as np
import os
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from netCDF4 import  Dataset
import sys
sys.path.append('../../shared/python/')

from util import make_folder
from helper.compute_monthly_stats import compute_monthly_stats
from helper.write_nc4 import write_sums_nc4, write_stats_nc4

import warnings; warnings.filterwarnings("ignore")

class obsfcstana_prep:
    
    def __init__(self, exp_list, start_time, end_time, obs_from=0):
        self.expdir_list = [item['expdir'] for item in exp_list]
        self.expid_list = [item['expid'] for item in exp_list]
        self.exptag_list = [item['exptag'] for item in exp_list]
        self.domain   = exp_list[0]['domain']
        self.start_time = start_time 
        self.end_time = end_time 
        self.var_list = ['obs_obs', 'obs_obsvar','obs_fcst','obs_fcstvar','obs_ana','obs_anavar']
        self.tilecoord = exp_list[0]['tilecoord']
        self.obsparam_list = [item['obsparam'] for item in exp_list]
        self.obs_from = obs_from

    def save_monthly_sum(self, outpath='./'):
        expdir_list  = self.expdir_list
        expid_list = self.expid_list
        exptag_list = self.exptag_list
        domain = self.domain
        var_list = self.var_list
        tc = self.tilecoord
        obsparam_list = self.obsparam_list
        start_time = self.start_time
        end_time = self.end_time

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
        
    def calculate_stats_fromsums(self, mo_path='./', write_to_nc=True, filename='./stats.nc4'):
        
        start_time = self.start_time
        end_time = self.end_time
        exptag_list = self.exptag_list
        
        # Variable list for computing sum and sum of squared
        var_list = self.var_list 

        # Read tilecoord and obsparam for tile and obs species information
        n_tile = self.tilecoord['N_tile']
        n_spec = len(self.obsparam_list[0])

        # Initialize statistical metrics 
        data_sum = {}
        data2_sum = {}
        N_data = np.zeros((n_tile, n_spec))
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
                    mN_data = nc.variables['N_data'][:]
                    moxf_sum = nc.variables['obsxfcst_sum'][:]
                    moxa_sum = nc.variables['obsxana_sum'][:]
                    mfxa_sum = nc.variables['fcstxana_sum'][:]
                    for var in var_list:
                        mdata_sum[var] = nc.variables[var+'_sum'][:]
                        mdata2_sum[var] = nc.variables[var+'2_sum'][:]
                       
                # Aggregate monthly data
                N_data += mN_data
                oxf_sum += moxf_sum
                oxa_sum += moxa_sum
                fxa_sum += mfxa_sum
               
                for var in var_list:
                    data_sum[var] += mdata_sum[var] 
                    data2_sum[var] += mdata2_sum[var]  
            else:
                raise FileNotFoundErroor(f"File {fout} does not exist")
                 
            current_month =current_month + relativedelta(months=1)

        # Compute the basic statistics after finishing time loop based on the accumulated data.
        data_mean ={}
        data2_mean = {}
        data_var = {}

        # calculation  
        for var in var_list:
            data_sum[var][N_data == 0] = np.nan
            data2_sum[var][N_data == 0] = np.nan
            
            data_mean[var]  = data_sum[var] / N_data
            data2_mean[var] = data2_sum[var] /N_data
            # var(x) = E[x2] - (E[x])^2
            data_var[var] = data2_mean[var] - data_mean[var]**2
            
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
        stats['N_data'] = N_data

        if write_to_nc:
            print('writing stats nc4 file: '+filename)
            write_stats_nc4(filename, stats)
            
        return stats


