import os
import numpy as np
import numpy.matlib as matlib
from scipy.io import loadmat
from netCDF4 import Dataset
from datetime import timedelta, datetime
from src.readwrite import read_GEOSIT, read_SCLF1C_nc, write_bin_SMOS_reg
from src.helper import gal_atm_correction, xy2hv
from src.helper.EASEv2 import EASEv2_latlon2ind, EASEv2_ind2latlon
import warnings; warnings.filterwarnings("ignore")

def preprocess_nc(smos_nc, outpath):

    # Code to preprocess SMOS data to the M36 EASEv2 grid.
    
    # Input:   1 SMOS-files: *.nc
    
    # Output:  1) M36, angular binned Tb w/ Faraday and geometric rotation
    #          2) M36, angular binned Tb w/ Faraday and geometric rotation,
    #             and removal of galactic, sun, moon,...
    #          3) M36, angular binned Tb w/ Faraday and geometric rotation,
    #             and removal of galactic, sun, moon,... (aux data)
    #             and removal of atmospheric correction (Peng et al., 2013)
    
    # *** SCLF1C ***
    
    # Qing Liu NASA/GSFC 2021 Convert matlab code to python	
    # Gabrielle De Lannoy - NASA/GSFC, based on code 21sept10

    #========================================================
    #  + Apply Tx/Ty conversion to Tbh/Tbv
    #  + Loop through SMOS pixels, assign EASE-indices,
    #    add (average) them per EASE-grid
    #    AND average them per degree bin
    #  + Sequential binary output (Fortran)
    
    #========================================================
    # GDL, 24Nov14:  - add additional output files (2, 3)
    # GDL, 21Aug15:  - fully remove any pointing to RWAPI-library
    #                - instead include xml-sparsing based on matlab routines
    #                - add preprocessor version in output
    #                - bugfix: 
    #                  + corrected spatial stdv of DGG values within EASE pixel  
    #                  + possibly corrected Tb_Sky radiation 
 
    overwrite=1
    
    # 1 = do overwrite to get the latest processed file at the same time stamp!!
    
    prep_version='v102'
    prep_version=int(prep_version[1:])
    
    # OUTPUT:
    #--------------------------------------------------------
    N_out_fields=16
    
    # Auxiliary info
    grid_name='M36'

    if grid_name == 'M36':
        N_ease_lat=406
        N_ease_lon=964
    else:
        raise RuntimeError('not ready for any other resolution than EASEv2 M36')
    
    write_ind_latlon='latlon'
    
    print('WRITE OUT...'+write_ind_latlon)

    # INPUT SPECS:
    N_sec_half_orbit=3239
    
    one_time_tag=14
    all_time_tag=2*one_time_tag + 1 + 1
    K_err=75
    
    max_K=375
    
    max_el=150000
    N_angle=200
    
    max_tot_el=max_el*N_angle
    out_min_res=30.0
    
    angle_step=1 
    
    start_angle=19.5
    
    end_angle=60.5
    
    inc_range=np.arange(start_angle,end_angle+0.1,angle_step)
    N_angle=np.size(inc_range) - 1
    inc_angle=np.arange(start_angle + angle_step / 2.0,end_angle - angle_step / 2.0 + 0.01,angle_step)
    if (N_angle != np.size(inc_angle)):
        raise RuntimeError('angle bins badly defined')
    
    correct_galaxy_only=0
    
    Asc='N'
    int_Asc=np.nan
    d_lon=np.full(max_tot_el,np.nan)
    d_lat=np.full(max_tot_el,np.nan)
    d_Tbh=np.full(max_tot_el,np.nan)
    d_Tbv=np.full(max_tot_el,np.nan)
    d_inc=np.full(max_tot_el,np.nan)
    d_Az=np.full(max_tot_el,np.nan)
    
    N_dat = N_ease_lat*N_ease_lon
    
    data_all=np.full((N_out_fields,N_dat,N_angle),0.)
    
    ind_i_tmp=matlib.repmat(np.array([np.arange(N_ease_lat)]).T,1,N_ease_lon)
    ind_j_tmp=matlib.repmat(np.array([np.arange(N_ease_lon)]),N_ease_lat,1)
    
    ind_i_all=np.reshape(ind_i_tmp,N_dat,order='F')
    ind_j_all=np.reshape(ind_j_tmp,N_dat,order='F')
    
    N_f=8
    data_gridded_bin=np.full((N_f,N_angle,N_ease_lat,N_ease_lon),0.)
    data_gridded_bin_sq=np.full((N_f,N_angle,N_ease_lat,N_ease_lon),0.)
    N_data_bin=np.full((N_f,N_angle,N_ease_lat,N_ease_lon),0.)
    
    # 2 datasets in the SCL product: dealt with in XY2HV
    print('READING')
       
    # read SMOS SCLF1C netcdf file
    TB_real,TB_imag,RA,Theta,Az,Fa,Ge,Snap_ID,t_smos_sec,lat,lon,flag_15,flag_16,mask_ok,Grid_ID,BT_count=read_SCLF1C_nc(smos_nc)
            
    fin =  Dataset(smos_nc,'r')
    Creator_Version=fin.getncattr('Fixed_Header:Source:Creator_Version')
    Asc=fin.getncattr('Variable_Header:Specific_Product_Header:Main_Info:Time_Info:Ascending_Flag')
    fname =  fin.getncattr('Fixed_Header:File_Name')
    idx = fname.index('SCLF1C_')
                       
    start_time = datetime(int(fname[idx+7:idx+11]),
                      int(fname[idx+11:idx+13]),
                      int(fname[idx+13:idx+15]),
                      int(fname[idx+16:idx+18]),
                      int(fname[idx+18:idx+20]),
                      int(fname[idx+20:idx+22]))
             
    end_time = datetime(int(fname[idx+23:idx+27]),
                      int(fname[idx+27:idx+29]),
                      int(fname[idx+29:idx+31]),
                      int(fname[idx+32:idx+34]),
                      int(fname[idx+34:idx+36]),
                      int(fname[idx+36:idx+38]))
        
    N_sec_date_diff =  end_time - start_time;
        
    if (N_sec_date_diff.total_seconds() < 0.9*N_sec_half_orbit):
        print('--This is an incomplete half orbit with '+ 
              str(N_sec_date_diff.total_seconds())+
              'sec instead of '+str(N_sec_half_orbit))
            
    t_mid = start_time + timedelta(seconds=np.round(N_sec_date_diff.total_seconds()/2))
    t_out = t_mid.replace(second=0)
        
    if np.round(t_mid.minute/out_min_res) < 1:
        t_out = t_out.replace(minute = 0)
    elif np.round(t_mid.minute/out_min_res) == 1:
        t_out = t_out.replace(minute = 30)
    else:
        t_out = t_out.replace(minute = 0)
        t_out = t_out + timedelta(seconds=3600);
          
    if int(Creator_Version) >= 344:
        alpha=Fa + Ge
    else:
        alpha=Ge - Fa
        
    del Fa,Ge
        
    d_lon,d_lat,d_inc,d_Az,end_id,d_Tbh,\
         d_Tbv,d_T3,d_T4,d_RAh,d_RAv,d_RA3,d_RA4,\
         d_Xorg,d_Yorg,d_RXYorg,d_IXYorg = xy2hv \
         (TB_real,TB_imag,RA,Theta,Az,alpha,Snap_ID,t_smos_sec,\
         lat,lon,flag_15,flag_16,mask_ok,Grid_ID,BT_count,nargout=17)
        
    d_Tbh[d_Xorg == 0]=np.nan
    d_Tbv[d_Yorg == 0]=np.nan
    d_T3[d_RXYorg == 0]=np.nan
    d_T4[d_IXYorg == 0]=np.nan
        
    dat_l=max(end_id)
        
    if Asc =='A':
        int_Asc=1
    else: 
        int_Asc=0
            
    N_points=len(d_Tbh)
    check=len(d_Tbv)
        
    if ((N_points != check) or (N_points != dat_l)):
        raise RuntimeError('error in variable dimensions coming out of rotation')
            
    print('Working on '+str(N_points)+' data points ; max='+str(max_tot_el))
    C=1

    # obtain GEOSIT latlon -> EASEv2 grid remapping info:
    # the mat file is pre-generated usign the "aux_grids_SMOS_prep" matlab script 

    m2d = loadmat(os.path.join('data','GEOSIT_to_EASEv2_M36.mat'))  
    vegcls_grid = m2d['vegcls_grid']
    NN_grid = m2d['NN_grid']
        
    T_air,T_surf,V_surf,P_surf = read_GEOSIT(t_out,vegcls_grid,NN_grid,d_lat, d_lon)
        
    Tb_ap_H,Tb_ap_V,Tb_BOA_H,Tb_BOA_V = \
                   gal_atm_correction(d_Tbh,d_Tbv,d_lat,d_lon,d_inc,d_Az,t_out,  \
                                      T_air,T_surf,P_surf,V_surf,C,nargout=4)

    if correct_galaxy_only:
        print('Galact corr (H):     maxdiff '+str(np.nanmax(Tb_ap_H - d_Tbh))+' mindiff '+str(np.nanmin(Tb_ap_H - d_Tbh)))
    else:
        print('Finished galact+atm correction')
        #print('Galact+atm corr (H): maxdiff '+str(np.nanmax(Tb_BOA_H - d_Tbh))+' mindiff '+str(np.nanmin(Tb_BOA_H - d_Tbh)))

    print('Project the SMOS data on a M36 grid...')
    ind_row,ind_col=EASEv2_latlon2ind(d_lat,d_lon,grid_name,nargout=2)
    if correct_galaxy_only:
        Tb_data=np.array([Tb_ap_H[np.arange(N_points)],
                     Tb_ap_V[np.arange(N_points)],
                     d_RAh[np.arange(N_points)],
                     d_RAv[np.arange(N_points)],
                     np.nan*d_T3[np.arange(N_points)],
                     np.nan*d_T4[np.arange(N_points)],
                     np.nan*d_RA3[np.arange(N_points)],
                     np.nan*d_RA4[np.arange(N_points)]])
    else:
        Tb_data=np.array([Tb_BOA_H[np.arange(N_points)],
                     Tb_BOA_V[np.arange(N_points)],
                     d_RAh[np.arange(N_points)],
                     d_RAv[np.arange(N_points)],
                     np.nan*d_T3[np.arange(N_points)],
                     np.nan*d_T4[np.arange(N_points)],
                     np.nan*d_RA3[np.arange(N_points)],
                     np.nan*d_RA4[np.arange(N_points)]])
        
    coord=np.array([ind_col[np.arange(N_points)],
               ind_row[np.arange(N_points)],
               d_inc[np.arange(N_points)]])
       
    for i in np.arange(N_points):
        ind_j=coord[0,i].astype('int')
        ind_i=coord[1,i].astype('int')
        ind_a=np.ceil((coord[2,i] - start_angle) / angle_step)
        ind_a=ind_a.astype('int')
            
        if np.any(np.logical_not(np.isnan(Tb_data[:,i]))):
            if ind_a > 0 and ind_a <= N_angle and \
                np.nanmax(Tb_data[:,i]) < max_K:
                ind_f=np.nonzero(np.logical_not(np.isnan(Tb_data[:,i])))[0]
                data_gridded_bin[ind_f,ind_a-1,ind_i,ind_j]=data_gridded_bin[ind_f,ind_a-1,ind_i,ind_j] + Tb_data[ind_f,i]
                data_gridded_bin_sq[ind_f,ind_a-1,ind_i,ind_j]=data_gridded_bin_sq[ind_f,ind_a-1,ind_i,ind_j] + Tb_data[ind_f,i]*Tb_data[ind_f,i]
                N_data_bin[ind_f,ind_a-1,ind_i,ind_j]=N_data_bin[ind_f,ind_a-1,ind_i,ind_j] + 1

    # Could refine the angle bins, filter out more RFI, i.e.
    # using K_err as maximum diff from mean,
    # and then reshape the data afterwards to get less bins
    N_data_bin[N_data_bin==0] = np.nan            
    data_gridded_bin=data_gridded_bin / N_data_bin
    data_gridded_bin_sq=data_gridded_bin_sq / N_data_bin
    data_gridded_bin_sq=np.sqrt((data_gridded_bin_sq - data_gridded_bin ** 2)*N_data_bin/ (N_data_bin - 1))

    # Set #points to zero, so that RFI-affected Tb 
    # won't show up in output later.
    N_data_bin[data_gridded_bin_sq > K_err]=0
    print('Finished the averaging\\n')
       
    for a in np.arange(N_angle):
        data_all[0,:,a]=np.reshape(data_gridded_bin[0,a,:,:],[N_dat],order='F')
        data_all[1,:,a]=np.reshape(data_gridded_bin[1,a,:,:],[N_dat],order='F')
        data_all[2,:,a]=np.reshape(data_gridded_bin_sq[0,a,:,:],[N_dat],order='F')
        data_all[3,:,a]=np.reshape(data_gridded_bin_sq[1,a,:,:],[N_dat],order='F')
        data_all[4,:,a]=np.reshape(N_data_bin[0,a,:,:],[N_dat],order='F')
        data_all[5,:,a]=np.reshape(N_data_bin[1,a,:,:],[N_dat],order='F')
        data_all[6,:,a]=np.reshape(data_gridded_bin[2,a,:,:],[N_dat],order='F')
        data_all[7,:,a]=np.reshape(data_gridded_bin[3,a,:,:],[N_dat],order='F')
        data_all[8,:,a]=np.reshape(data_gridded_bin[4,a,:,:],[N_dat],order='F')
        data_all[9,:,a]=np.reshape(data_gridded_bin[5,a,:,:],[N_dat],order='F')
        data_all[10,:,a]=np.reshape(data_gridded_bin_sq[4,a,:,:],[N_dat],order='F')
        data_all[11,:,a]=np.reshape(data_gridded_bin_sq[5,a,:,:],[N_dat],order='F')
        data_all[12,:,a]=np.reshape(N_data_bin[4,a,:,:],[N_dat],order='F')
        data_all[13,:,a]=np.reshape(N_data_bin[5,a,:,:],[N_dat],order='F')
        data_all[14,:,a]=np.reshape(data_gridded_bin[6,a,:,:],[N_dat],order='F')
        data_all[15,:,a]=np.reshape(data_gridded_bin[7,a,:,:],[N_dat],order='F')
          
        ind=np.union1d(np.nonzero(data_all[4,:,a] > 0)[0],np.nonzero(data_all[5,:,a] > 0)[0])
    
    data_all[np.isnan(data_all)]= -999.0
    stamp1=t_out.strftime("%Y")
    stamp2=t_out.strftime("%m")
    stamp3=t_out.strftime("%d")
    stamp4=t_out.strftime("%H")
    stamp5=t_out.strftime("%M")
         
    if not os.path.exists(outpath+'/'+stamp1+stamp2):
        os.makedirs(outpath+'/'+stamp1+stamp2)
            
    if correct_galaxy_only:
        out_filename=outpath+'/'+stamp1+stamp2+'/SMOS_reg_nosky_Tb_'+stamp1+stamp2+stamp3+'_'+stamp4+stamp5+'_'+Asc+'.bin'
    else:
        out_filename=outpath+'/'+stamp1+stamp2+'/SMOS_reg_Tb_'+stamp1+stamp2+stamp3+'_'+stamp4+stamp5+'_'+Asc+'.bin'
    
    if ind.size > 0:
        if write_ind_latlon == 'ind':
            write_bin_SMOS_reg(out_filename,ind_j_all[ind],ind_i_all[ind],inc_angle,data_all[:,ind,:],int_Asc,Creator_Version,prep_version,start_time,end_time,overwrite,N_out_fields,write_ind_latlon,'SCLF1C')
        else:
            if write_ind_latlon == 'latlon':
                lat_out,lon_out=EASEv2_ind2latlon(ind_i_all[ind],ind_j_all[ind],grid_name,nargout=2)
                write_bin_SMOS_reg(out_filename,lon_out,lat_out,inc_angle,data_all[:,ind,:],int_Asc,int(Creator_Version),prep_version,start_time,end_time,overwrite,N_out_fields,write_ind_latlon,'SCLF1C')
            else:
                raise RuntimeError('Unknown format of indices / latlon output')
    else:
        print('Zero valid data, NO output file written')   

