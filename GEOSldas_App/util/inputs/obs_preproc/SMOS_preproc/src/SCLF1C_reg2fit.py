from datetime import datetime,timedelta
import numpy as np
import os
import sys
from src.readwrite import read_bin_SMOS_reg, write_bin_SMOS_reg

def SCLF1C_reg2fit(reg_path,start_time,end_time,Asc):
    #Read preprocessed 'regular' SMOS Tb data - all angles
    #Interpolation / fit -> extract 40^o angle, and all other angles
    #========================================================
    
    in_path=reg_path
    out_path=reg_path.replace('_reg_','_fit_')
    #file_in_prefix  = 'SMOS_reg_nosky_noatm_Tb_';
    #file_out_prefix = 'SMOS_fit_nosky_noatm_Tb_';
    file_in_prefix='SMOS_reg_Tb_'
    file_out_prefix='SMOS_fit_Tb_'
 
    # half a hour
    dtstep=1800 
    #----------------------------------
    TbHVconstrained=''
    #TbHVconstrained = '_TbHVconstrained'; #constrained fitting of multiple Tb-variables, or empty if individual
        
    #EXPONENTIAL fit: y = b0 + b1.exp(-b3.x)
    #tunable coefficients appear non-linearly -> non-linear fitting takes iteration time
    #----------------------------------
    #fit_type = 'exp1';
        
    #QUADRATIC fit: y =ax^2 + bx + c
    #tunable coefficients appear linearly -> quick matrix inversion
    #----------------------------------
    fit_type='poly2'
        
    read_ind_latlon='latlon'
    write_ind_latlon='latlon'
        
    overwrite=1
    
    #------------------------------------------------------
    # QC
    #------------------------------------------------------
        
    N_gt30_lt50=10     #number of data at angles>30 and <50^o
        
    N_points_tot=15    #total number of angles used in the fit
        
    min_N_ang=5        #request at least 5 angles available at either side
                       #of each interpolated value
       
    #remove data that indicate of RFI:
    cutoff_weight=2 / (7 ** 2)  #eliminate data with excessive Tb variability
        
    cutoff_K=320
        
    #------------------------------------------------------
        
    N_out_fields=16
    field_tags=['Tb_H','Tb_V','std_h','std_v','c_h','c_v','RAh','RAv',\
                'Tb3','Tb4','std_3','std_4','c_3','c_4','RA3','RA4']
    
    #Variables to be fitted
    #-------------------------
    #Tb_ind  = [0 1  8 9];
    #We really do not know the shape of T3, T4 - do not try to fit for now
        
    Tb_ind=[0,1] #Tbh, Tbv
        
    # std_ind = Tb_ind + 2;
    # c_ind   = Tb_ind + 4;
    # Ra_ind  = Tb_ind + 6;
        
    data_product='SCLF1C'
        
    date_time_new=start_time
    
    no_data_value=-999.0
    no_data_tol=0.0001
        
    # Go through loop of orbits
    while date_time_new < end_time:
        # augment date_time and t_ind
        date_time_old=date_time_new
        date_time_new=date_time_old + timedelta(seconds=dtstep)
    
        fname = in_path+date_time_new.strftime('%Y%m')+'/'+ \
                file_in_prefix+date_time_new.strftime('%Y%m%d_%H%M')+\
                Asc+'.bin'
            
        if not os.path.exists(out_path+'/SMOS_fit_'+fit_type+TbHVconstrained+'/'+ \
                date_time_new.strftime('%Y%m')):
            os.makedirs(out_path+'/SMOS_fit_'+fit_type+TbHVconstrained+'/'+ \
                date_time_new.strftime('%Y%m'))  
                
        if os.path.isfile(fname):
            out_filename=out_path+'/SMOS_fit_'+fit_type+TbHVconstrained+'/'+ \
                date_time_new.strftime('%Y%m')+'/'+file_out_prefix+ \
                date_time_new.strftime('%Y%m%d_%H%M')+Asc+'.bin'
    
            time=date_time_new.strftime('%Y%m%d%H%M')       
            
            data,inc_ang,col_ind,row_ind,asc_flag,version,prep_version,o_start_time,o_end_time,N_grid=\
            read_bin_SMOS_reg(fname,N_out_fields,read_ind_latlon,data_product,nargout=10)
                
            #inc_ang=inc_ang.T
            data_out=data.copy()

            if len(data_out.shape) < 3: 
                continue
           
            for v in range(len(Tb_ind)):
                data_out[Tb_ind[v],:,:]=np.nan
                data_out[Tb_ind[v] + 2,:,:]=np.nan
                data_out[Tb_ind[v] + 4,:,:]=np.nan
                data_out[Tb_ind[v] + 6,:,:]=np.nan
            #--------------------------------------------------------------------
            #weighted interpolation
            #weight = number of points in average / stdv
            #--------------------------------------------------------------------
            data[np.abs(data - no_data_value) < no_data_tol]=np.nan
            
            #for all grid cells, overwrite the original 'data'    
            #with new 'data' containing fitted information
            for i in range(len(col_ind)):
                if (~np.isnan(data[:,i,:])).any():
                    if TbHVconstrained == '':
                        #for each of the 4 or 2 Tb-variables individually
                        for v in range(len(Tb_ind)):
                            #1) actual Tb
                            Tb_data=data[Tb_ind[v],i,:]
                            Tb_data[Tb_data < 0]=np.nan
                                
                            #stdv+RA
                            #total_prior_uncert = sqrt(squeeze(data(Tb_ind(v)+2,i,:)).^2 + squeeze(data(Tb_ind(v)+6,i,:)).^2);
                            total_prior_uncert=data[Tb_ind[v]+2,i,:]**2 + \
                                               data[Tb_ind[v]+6,i,:]**2
                            #Each angle has a different number of contributing DGG cells, 
                            #because it is collecting data over a [x-0.5 x+0.5]-window and 
                            #from different instants (seconds) in time (fwd, backward looking)
                            weights=np.sqrt(data[Tb_ind[v]+4,i,:] / total_prior_uncert)
                            good=np.nonzero(~np.isnan(Tb_data) * (Tb_data < cutoff_K) * \
                                            (weights>cutoff_weight) * ~np.isnan(weights) )[0]
                                
                            Tb_data=Tb_data[good]
                            weights=weights[good]
                            x=inc_ang[good]
                            
                            if ((x[(x >= 30) * (x<= 50)].size >= N_gt30_lt50 )* 
                                (x.size >= N_points_tot) * \
                                ((fit_type == 'exp1') + \
                                 ((fit_type == 'poly2') * \
                                  (x[x <= 40].size >= min_N_ang) * \
                                  (x[x >= 40].size >= min_N_ang)))):
                                    #quadratic fit
                                    if fit_type == 'poly2':
                                        #Quadratic fit, no constraints.
                                        p_w=np.polyfit(x,Tb_data,2, w=weights)
                                        Tb_data_fit=np.polyval(p_w,inc_ang)
                                        #available at either side of the interpolation
                                        count_angle=np.cumsum(1 + 0.*good)
                                        good_4=np.nonzero(np.logical_and(count_angle > min_N_ang, (count_angle + min_N_ang) <= count_angle[-1]))[0]
                                        only_good=good[good_4]
                                    else:
                                        sys.exit('fitfunction not available')
                                    #==> replace Radiometric Accuracy w/ sqrt(average(Ra)^2)
                                    RA_data=data[Tb_ind[v]+6,i,:]
                                    RA_data[RA_data < 0]=np.nan
                                    RA_data=RA_data[only_good]
                                    RA_data_new=np.sqrt(np.mean(RA_data ** 2))
                                    stdv_data_new=np.sqrt(np.mean((Tb_data[good_4] - Tb_data_fit[only_good]) ** 2))
                                    #counts and stdv (equal for all angles)
                                    data_out[Tb_ind[v],i,:]=np.nan
                                    data_out[Tb_ind[v],i,only_good]=Tb_data_fit[only_good]
                                    data_out[Tb_ind[v] + 2,i,:]=stdv_data_new
                                    data_out[Tb_ind[v] + 4,i,:]=only_good.size
                                    data_out[Tb_ind[v] + 6,i,:]=RA_data_new
                            else:
                                data_out[Tb_ind[v],i,:]=np.nan
                                data_out[Tb_ind[v] + 2,i,:]=np.nan
                                data_out[Tb_ind[v] + 4,i,:]=np.nan
                                data_out[Tb_ind[v] + 6,i,:]=np.nan
                    else:
                        sys.exit('Revise fitting options')
                        
            data=data_out
            if (~np.isnan(data)).any():
                data[np.isnan(data)]=no_data_value 
                write_bin_SMOS_reg(out_filename,col_ind,row_ind,inc_ang,data,\
                                      asc_flag,version,prep_version,o_start_time,\
                                          o_end_time,overwrite,N_out_fields,\
                                              write_ind_latlon,'SCLF1C')


