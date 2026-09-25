import numpy as np
from numpy import matlib 

def invMR4L(alpha = None,*args,**kwargs):
 
    cosdalpha = np.cos(np.deg2rad(alpha))
    sindalpha = np.sin(np.deg2rad(alpha))
    cosd2alpha = np.cos(np.deg2rad(2*alpha))
    sind2alpha = np.sin(np.deg2rad(2*alpha))
    
    MR4 = np.full([len(alpha),4,4],np.nan)
    
    MR4[:,0,0] = cosdalpha ** 2
    MR4[:,0,1] = sindalpha ** 2
    MR4[:,0,2] = -cosdalpha * sindalpha
    MR4[:,0,3] = 0
    MR4[:,1,0] = sindalpha ** 2
    MR4[:,1,1] = cosdalpha ** 2
    MR4[:,1,2] = cosdalpha*sindalpha
    MR4[:,1,3] = 0
    MR4[:,2,0] = sind2alpha
    MR4[:,2,1] = -sind2alpha
    MR4[:,2,2] = cosd2alpha
    MR4[:,2,3] = 0
    MR4[:,3,0] = 0
    MR4[:,3,1] = 0
    MR4[:,3,2] = 0
    MR4[:,3,3] = 1
    
    IMR4 = np.full(MR4.shape,np.nan)
    
    IMR4[:,0,0]=MR4[:,1,2]*MR4[:,2,3]*MR4[:,3,1] - MR4[:,1,3]*MR4[:,2,2]*MR4[:,3,1] + MR4[:,1,3]*MR4[:,2,1]*MR4[:,3,2] - MR4[:,1,1]*MR4[:,2,3]*MR4[:,3,2] - MR4[:,1,2]*MR4[:,2,1]*MR4[:,3,3] + MR4[:,1,1]*MR4[:,2,2]*MR4[:,3,3]
    IMR4[:,0,1]=MR4[:,0,3]*MR4[:,2,2]*MR4[:,3,1] - MR4[:,0,2]*MR4[:,2,3]*MR4[:,3,1] - MR4[:,0,3]*MR4[:,2,1]*MR4[:,3,2] + MR4[:,0,1]*MR4[:,2,3]*MR4[:,3,2] + MR4[:,0,2]*MR4[:,2,1]*MR4[:,3,3] - MR4[:,0,1]*MR4[:,2,2]*MR4[:,3,3]
    IMR4[:,0,2]=MR4[:,0,2]*MR4[:,1,3]*MR4[:,3,1] - MR4[:,0,3]*MR4[:,1,2]*MR4[:,3,1] + MR4[:,0,3]*MR4[:,1,1]*MR4[:,3,2] - MR4[:,0,1]*MR4[:,1,3]*MR4[:,3,2] - MR4[:,0,2]*MR4[:,1,1]*MR4[:,3,3] + MR4[:,0,1]*MR4[:,1,2]*MR4[:,3,3]
    IMR4[:,0,3]=MR4[:,0,3]*MR4[:,1,2]*MR4[:,2,1] - MR4[:,0,2]*MR4[:,1,3]*MR4[:,2,1] - MR4[:,0,3]*MR4[:,1,1]*MR4[:,2,2] + MR4[:,0,1]*MR4[:,1,3]*MR4[:,2,2] + MR4[:,0,2]*MR4[:,1,1]*MR4[:,2,3] - MR4[:,0,1]*MR4[:,1,2]*MR4[:,2,3]
    IMR4[:,1,0]=MR4[:,1,3]*MR4[:,2,2]*MR4[:,3,0] - MR4[:,1,2]*MR4[:,2,3]*MR4[:,3,0] - MR4[:,1,3]*MR4[:,2,0]*MR4[:,3,2] + MR4[:,1,0]*MR4[:,2,3]*MR4[:,3,2] + MR4[:,1,2]*MR4[:,2,0]*MR4[:,3,3] - MR4[:,1,0]*MR4[:,2,2]*MR4[:,3,3]
    IMR4[:,1,1]=MR4[:,0,2]*MR4[:,2,3]*MR4[:,3,0] - MR4[:,0,3]*MR4[:,2,2]*MR4[:,3,0] + MR4[:,0,3]*MR4[:,2,0]*MR4[:,3,2] - MR4[:,0,0]*MR4[:,2,3]*MR4[:,3,2] - MR4[:,0,2]*MR4[:,2,0]*MR4[:,3,3] + MR4[:,0,0]*MR4[:,2,2]*MR4[:,3,3]
    IMR4[:,1,2]=MR4[:,0,3]*MR4[:,1,2]*MR4[:,3,0] - MR4[:,0,2]*MR4[:,1,3]*MR4[:,3,0] - MR4[:,0,3]*MR4[:,1,0]*MR4[:,3,2] + MR4[:,0,0]*MR4[:,1,3]*MR4[:,3,2] + MR4[:,0,2]*MR4[:,1,0]*MR4[:,3,3] - MR4[:,0,0]*MR4[:,1,2]*MR4[:,3,3]
    IMR4[:,1,3]=MR4[:,0,2]*MR4[:,1,3]*MR4[:,2,0] - MR4[:,0,3]*MR4[:,1,2]*MR4[:,2,0] + MR4[:,0,3]*MR4[:,1,0]*MR4[:,2,2] - MR4[:,0,0]*MR4[:,1,3]*MR4[:,2,2] - MR4[:,0,2]*MR4[:,1,0]*MR4[:,2,3] + MR4[:,0,0]*MR4[:,1,2]*MR4[:,2,3]
    IMR4[:,2,0]=MR4[:,1,1]*MR4[:,2,3]*MR4[:,3,0] - MR4[:,1,3]*MR4[:,2,1]*MR4[:,3,0] + MR4[:,1,3]*MR4[:,2,0]*MR4[:,3,1] - MR4[:,1,0]*MR4[:,2,3]*MR4[:,3,1] - MR4[:,1,1]*MR4[:,2,0]*MR4[:,3,3] + MR4[:,1,0]*MR4[:,2,1]*MR4[:,3,3]
    IMR4[:,2,1]=MR4[:,0,3]*MR4[:,2,1]*MR4[:,3,0] - MR4[:,0,1]*MR4[:,2,3]*MR4[:,3,0] - MR4[:,0,3]*MR4[:,2,0]*MR4[:,3,1] + MR4[:,0,0]*MR4[:,2,3]*MR4[:,3,1] + MR4[:,0,1]*MR4[:,2,0]*MR4[:,3,3] - MR4[:,0,0]*MR4[:,2,1]*MR4[:,3,3]
    IMR4[:,2,2]=MR4[:,0,1]*MR4[:,1,3]*MR4[:,3,0] - MR4[:,0,3]*MR4[:,1,1]*MR4[:,3,0] + MR4[:,0,3]*MR4[:,1,0]*MR4[:,3,1] - MR4[:,0,0]*MR4[:,1,3]*MR4[:,3,1] - MR4[:,0,1]*MR4[:,1,0]*MR4[:,3,3] + MR4[:,0,0]*MR4[:,1,1]*MR4[:,3,3]
    IMR4[:,2,3]=MR4[:,0,3]*MR4[:,1,1]*MR4[:,2,0] - MR4[:,0,1]*MR4[:,1,3]*MR4[:,2,0] - MR4[:,0,3]*MR4[:,1,0]*MR4[:,2,1] + MR4[:,0,0]*MR4[:,1,3]*MR4[:,2,1] + MR4[:,0,1]*MR4[:,1,0]*MR4[:,2,3] - MR4[:,0,0]*MR4[:,1,1]*MR4[:,2,3]
    IMR4[:,3,0]=MR4[:,1,2]*MR4[:,2,1]*MR4[:,3,0] - MR4[:,1,1]*MR4[:,2,2]*MR4[:,3,0] - MR4[:,1,2]*MR4[:,2,0]*MR4[:,3,1] + MR4[:,1,0]*MR4[:,2,2]*MR4[:,3,1] + MR4[:,1,1]*MR4[:,2,0]*MR4[:,3,2] - MR4[:,1,0]*MR4[:,2,1]*MR4[:,3,2]
    IMR4[:,3,1]=MR4[:,0,1]*MR4[:,2,2]*MR4[:,3,0] - MR4[:,0,2]*MR4[:,2,1]*MR4[:,3,0] + MR4[:,0,2]*MR4[:,2,0]*MR4[:,3,1] - MR4[:,0,0]*MR4[:,2,2]*MR4[:,3,1] - MR4[:,0,1]*MR4[:,2,0]*MR4[:,3,2] + MR4[:,0,0]*MR4[:,2,1]*MR4[:,3,2]
    IMR4[:,3,2]=MR4[:,0,2]*MR4[:,1,1]*MR4[:,3,0] - MR4[:,0,1]*MR4[:,1,2]*MR4[:,3,0] - MR4[:,0,2]*MR4[:,1,0]*MR4[:,3,1] + MR4[:,0,0]*MR4[:,1,2]*MR4[:,3,1] + MR4[:,0,1]*MR4[:,1,0]*MR4[:,3,2] - MR4[:,0,0]*MR4[:,1,1]*MR4[:,3,2]
    IMR4[:,3,3]=MR4[:,0,1]*MR4[:,1,2]*MR4[:,2,0] - MR4[:,0,2]*MR4[:,1,1]*MR4[:,2,0] + MR4[:,0,2]*MR4[:,1,0]*MR4[:,2,1] - MR4[:,0,0]*MR4[:,1,2]*MR4[:,2,1] - MR4[:,0,1]*MR4[:,1,0]*MR4[:,2,2] + MR4[:,0,0]*MR4[:,1,1]*MR4[:,2,2]
    
    return IMR4

def xy2hv(a_BT_Real=None,a_BT_Imag=None,a_RA=None,a_Theta=None,\
                           a_Az=None,a_alpha=None, a_Snap_ID=None,a_t_smos_sec=None,\
                               lat=None,lon=None,a_flag_15=None,a_flag_16=None,\
                                   mask_ok=None,dgg_list=None,BT_count=None,\
                                       in_dgg=None, *args, **kwargs):
    '''
    # This function do the rotation of SMOS data from the antenna reference
    # (XY) to Earth reference (HV) - Top Of Atmosphere.
    # 
    # INPUTS :
    #          = output from call to TB nc file
        
    # OUTPUT : 
    #          = specific variables in a long 1D array
        
    # Interpolation and rotation can only be done if the number of available
    # angles is greater than 6 for full pol. If the roation
    # cannot be performed, the strucutre will be filled with NaN and a message 
    # will be displayed on screen.
    # 
    # Based on:
    # ==========
    # Authors  : Delphine Leroux (delphine.leroux@cesbio.cnes.fr)
    # Creation : 29/10/2010 (v1.1)
    # Version  : 1.8
    # Comments : - this fnunction has been tested on 64-bit Linux with RWAPI 
    # v1.3 and Matlab R2010a
    #            - in dual polarization, the user should be aware that the
    #            rotation matrix MR2 is not invertible around +-45° so the
    #            inversion might not succeed or will give very high values of
    #            brightness temperature. Therefore, the user can set a
    #            threshold or can filter the results.
    # v1.1 : no need to put the polarization as an input (directly read in the
    # filename) - 29/10/2010
    # v1.2 : in the full pol part, the beginning of the interpolation has been
    # changed (only the time field gives the order to sort all the other
    # fields) & lat, lon, lat and mask are stored - 02/11/2010
    # v1.3 : a reciprocal condition number threshold has been added such that 
    # if the rcond number of MR2 is less than this threshold, the MR2 matrix is
    # not inverted - only for DUAL pol - 02/11/2010
    # v1.4 : The flag for each measurement is now stored and can be used to
    # look at the RFI flag & a test has been added to check if the DGG of
    # interest are located in the input product - 04/11/2010
    # v1.5 : Bug fixed (lat/long were not stored correctly) & problem with the
    # alpha angle used to do the rotation fixed - 08/11/2010
    # v1.6 : minor changes to gain some speed - 10/11/2010
    # V1.7 : convention of the rotation angle changing depending on the version
    # of L1OP. For v344 or later a=Fa+Ge and before v344 a=Fa-Ge - 15/11/2010
    # v1.8 : The azimuth angles are now stored and the flag are stored
    # differently (now it is a matrix) - 17/11/2010
    # CESBIO - 2010
        
    # UPDATED:
    #==========
    # Gabrielle De Lannoy, NASA/GMAO, 17dec10: pre-allocated memory,
    #                                          removed loops,         
    #                                          updated for F pol only
    # Gabrielle De Lannoy, NASA/GMAO, 08mar11: take mex-output,
    #                                          insert optimization (FC)
    # Gabrielle De Lannoy, NASA/GMAO, 24mar11: reduce data coming in from reader
    #                                          by a priori selection on flags,
    #                                          limit sort/find/loop - commands,
    #                                          correct RA (to match v1.10)
    # Gabrielle De Lannoy, NASA/GMAO, 02dec14: added azimuth output
    # Gabrielle De Lannoy, NASA/GMAO, 24aug15: added check on interpolation
    #=====================================================
    '''
    #nargin = XY2HV_F_mex_opt_2dec14.nargin
    
    # 10 K = 2*4K +2K
    max_dev=10
    
    # don't smooth when there is a gap of more than 1 degree
    theta_threshold=1
    
    #if nargin < 19 :
    in_dgg = dgg_list
    ia = np.arange(len(in_dgg))
    
    if in_dgg.size >= 1:
        #######################        
        ## FULL POLARIZATION ##        
        #######################
        # Initializing the output structure
        
        # number of *valid* (no RFI, ALIAS-free) snapshots/inc angles per grid cell
        max_el = 250
        max_tot_el = len(in_dgg) * max_el
        
        out_lon= np.full([max_tot_el], np.nan)  
        out_lat= np.full([max_tot_el], np.nan)
        out_inc= np.full([max_tot_el], np.nan)
        out_Az= np.full([max_tot_el], np.nan)  
        out_Tbh= np.full([max_tot_el], np.nan)
        out_Tbv= np.full([max_tot_el], np.nan)
        out_T3= np.full([max_tot_el], np.nan)  
        out_T4= np.full([max_tot_el], np.nan)
        out_RAh= np.full([max_tot_el], np.nan)
        out_RAv= np.full([max_tot_el], np.nan)  
        out_RA3= np.full([max_tot_el], np.nan)
        out_RA4= np.full([max_tot_el], np.nan)
        
        out_Xorg= np.full([max_tot_el], 0.)  
        out_Yorg= np.full([max_tot_el], 0.)
        out_RXYorg= np.full([max_tot_el], 0.)
        out_IXYorg= np.full([max_tot_el], 0.)

        end_id_list = np.full([len(in_dgg)], np.nan)
        
        start_id = 0
        
        for idx_in_dgg in np.arange(len(in_dgg)):
            
            idx_dgg = ia[idx_in_dgg]
            TB_count = BT_count[idx_dgg]
            
            # only select land points,
            # and limit data to those points at a distance more than 40km away
            # from 'coastline' (i.e. water) 
            # => taken together in mask_ok (alias and RFI flags are
            # accounted for in reader)
            
            if (TB_count > 6 and mask_ok[idx_dgg] and lat[idx_dgg] > -60):

                flag_15 = a_flag_15[idx_dgg, :TB_count]
                flag_16 = a_flag_16[idx_dgg, :TB_count]
                BT_Real = a_BT_Real[idx_dgg, :TB_count]
                BT_Imag = a_BT_Imag[idx_dgg, :TB_count]
                
                RA      = a_RA[idx_dgg,:TB_count]
                Theta   = a_Theta[idx_dgg,:TB_count]
                Az      = a_Az[idx_dgg,:TB_count]
                alpha   = a_alpha[idx_dgg,:TB_count]
                Snap_ID = a_Snap_ID[idx_dgg,:TB_count]
                t_smos_sec = a_t_smos_sec[idx_dgg,:TB_count]
                
                # Only consider unique incidence angle, unique acquisition time and unique
                # snapshot ID even for mixed snapshot (XX and XY or YY and XY)
                # Unique and in the same order as found (sort(index))
                dum,b_tsmos = np.unique(t_smos_sec,return_index=True)
                
                b_tsmos = b_tsmos[np.logical_not(np.isnan(dum))]
                ind = np.sort(b_tsmos)
                t_smos_sec_uniq = t_smos_sec[ind]
                Theta_uniq = Theta[ind]
                Az_uniq = Az[ind]
                alpha_uniq = alpha[ind]
                Snap_ID_uniq = Snap_ID[ind]
                
                if Snap_ID_uniq.size == 0:
                    
                    end_id = start_id + 1
                    end_id_list[idx_in_dgg] = end_id 
                    out_lon[start_id:end_id] = np.nan
                    out_lat[start_id:end_id] = np.nan
                    out_Tbh[start_id:end_id] = np.nan
                    out_Tbv[start_id:end_id] = np.nan
                    out_T3[start_id:end_id] = np.nan
                    out_T4[start_id:end_id] = np.nan
                    out_RAh[start_id:end_id] = np.nan
                    out_RAv[start_id:end_id] = np.nan
                    out_RA3[start_id:end_id] = np.nan
                    out_RA4[start_id:end_id] = np.nan
                    
                    out_Xorg[start_id:end_id] = np.nan
                    out_Yorg[start_id:end_id] =  np.nan
                    out_RXYorg[start_id:end_id] = np.nan
                    out_IXYorg[start_id:end_id] = np.nan
                    out_inc[start_id:end_id] = np.nan
                    out_Az[start_id:end_id] = np.nan
                    start_id = end_id 
                    
                else:
                    
                    # np.nonzero return a tuple, thus need to specify [0] for 1-D idx
                    idx_TBxx = np.nonzero(np.logical_and(flag_15==0, flag_16 == 0))[0]
                    idx_TByy = np.nonzero(np.logical_and(flag_15==0, flag_16 == 1))[0]
                    idx_TBxy = np.nonzero(np.logical_or(np.logical_and(flag_15==1,flag_16==0), \
                                     np.logical_and(flag_15==1,flag_16==1)))[0]
                    
                    # TBxy and Flags_final = output; since we filter out for
                    # the flags in the step above, there is no need to have flags in
                    # output
                        
                    # TBxy : 1-TBxx 2-TByy 3-Re(TBxy) 4-Im(TBxy) 5-Theta 6-SnapID 7-t_smos_sec
                    #        8-RA_TBxx  9-RA_TByy  10-RA_TBxy
                 
                    TBxy = np.full([len(Snap_ID_uniq),10],np.nan)
                    
                    TBxy[:,4] = Theta_uniq
                    TBxy[:,5] = Snap_ID_uniq
                    TBxy[:,6] = t_smos_sec_uniq
                    
                    for i in np.arange(len(idx_TBxx)):
                        idx_snap = np.nonzero(Snap_ID_uniq == Snap_ID[idx_TBxx[i]])[0]
                        TBxy[idx_snap,0] = BT_Real[idx_TBxx[i]]
                        TBxy[idx_snap,7] = RA[idx_TBxx[i]]
                    for i in np.arange(len(idx_TByy)):
                        idx_snap = np.nonzero(Snap_ID_uniq == Snap_ID[idx_TByy[i]])[0]
                        TBxy[idx_snap,1] = BT_Real[idx_TByy[i]]
                        TBxy[idx_snap,8] = RA[idx_TByy[i]]
                    for i in np.arange(len(idx_TBxy)):
                        idx_snap = np.nonzero(Snap_ID_uniq == Snap_ID[idx_TBxy[i]])[0]
                        TBxy[idx_snap,2] = BT_Real[idx_TBxy[i]]
                        TBxy[idx_snap,3] = BT_Imag[idx_TBxy[i]]
                        TBxy[idx_snap,9] = RA[idx_TBxy[i]]
                        
                    #------------------------------------------------------------------
                    # Eliminate outlier data before interpolation
                    #------------------------------------------------------------------
                    # For that we need to do moving average per angle, because the TB
                    # will not necessarily follow a smooth path over time, and there
                    # may be unforeseen angle-gaps when smoothing over snapshots...
                    # It would be harder to select sub-windows of 'similar' snapshots;
                    # now we smooth only over continuous sets of angles
                            # ==> this could be easily reverted by limiting the
                    # snapshot-windows to periods with minimal time differences...
                        
                    TBxy_thetasorted = TBxy[np.argsort(TBxy[:,4]),:]
                    theta_1 = TBxy_thetasorted[0:-1,4]
                    theta_2 = TBxy_thetasorted[1:,4]
                    smooth_end = np.nonzero(np.abs(theta_1-theta_2) > theta_threshold)[0]
                    
                    smooth_end = np.concatenate(([-1],smooth_end,[len(theta_1)+1]),axis=None)
                                    
                    #Per continuous set of angles, do the smoothing and throw out outliers
                    N = 5
                    
                    for i in np.arange(len(smooth_end)-1):
                        
                        # Nothing is done on outliers 
                        # if there is no more than N consecutive angles
                        # if np.abs((smooth_end[i]+1) - smooth_end[i+1]) >= 5:
                        if np.abs((smooth_end[i]+1) - smooth_end[i+1]) > 5:
                           
                            tmp_TBxy = TBxy_thetasorted[smooth_end[i]+1:smooth_end[i+1],:4]
                            tmp = np.nancumsum(tmp_TBxy,0)
                            count_tmp = np.cumsum(np.logical_not(np.isnan(tmp_TBxy)),0)
                            
                            #Initial rows of tmp may contain 0, 
                            #if there happened to be leading nan-values for some fields.
                            #You would hope to have at least 1 good value in the first
                            #window of N angles - if not, look for the first non-zero
                            #element and replace zero cumsum by this first value
                            #assign window-mean to the middle element
                            #if count_tmp(N+1:end,:) - count_tmp(1:end-N,:) happens to be
                            #zero, then NaN results in the window mean, which is fine,
                            #because there won't be any data to check in this window anyway
                            #In case there was no single good value in the first N angles:
                            #Replace leading zeros (~nan) with first non-zero element
                            
                            for jj in np.arange(4):
                                if (count_tmp[0,jj] == 0):
                                    # index of first non-zero only
                                    ind=np.nonzero(tmp[:,jj] != 0)[0]
                                    if ind.size > 0:
                                        tmp[np.arange(ind[0]),jj]=tmp[ind[0],jj]
                                        count_tmp[np.arange(ind[0]),jj]=1
                                    
                            aux = tmp[N,:] / count_tmp[N,:] #for first window middle
                            tmp[int(np.ceil(N/2)):-int(np.floor(N/2)),:] = \
                                (tmp[N:,:]-tmp[0:-N,:]) / (count_tmp[N:,:]-count_tmp[0:-N,:])
                            tmp[int(np.ceil(N/2))-1,:]=aux
                            
                            # Assign constant (mean) values to tails of the moving window
                            tmp[0:int(np.ceil(N/2))-1,:]=np.matlib.repmat(tmp[int(np.ceil(N/2))-1,:], int(np.ceil(N/2))-1,1)
                            tmp[-int(np.floor(N/2)):,:]=np.matlib.repmat(tmp[-int(np.floor(N/2))-1,:],int(np.floor(N/2)),1)
                            
                            # Remove outliers
                            tmp_TBxy[np.abs(tmp_TBxy-tmp) > max_dev] = np.nan
                            TBxy_thetasorted[smooth_end[i]+1:smooth_end[i+1],0:4]=tmp_TBxy
                            
                            del tmp, count_tmp
                    
                    # sort back on time
                    TBxy = TBxy_thetasorted[np.argsort(TBxy_thetasorted[:,6]),:]
                    
                    #--------------------------------------------------------------------
    
                    # Locate the X, Y and XY polarizations
                    bool_Flags_X = np.logical_not(np.isnan(TBxy[:,0]))
                    bool_Flags_Y = np.logical_not(np.isnan(TBxy[:,1]))
                    bool_Flags_XY = np.logical_and(np.logical_not(np.isnan(TBxy[:,2])),\
                                                   np.logical_not(np.isnan(TBxy[:,3])))
                    
                    # To keep track of where original data were before interpolation   
                    bool_Flags_RXY = np.logical_not(np.isnan(TBxy[:,2]))
                    bool_Flags_IXY = np.logical_not(np.isnan(TBxy[:,3]))
                    
                    # Identify if the possibleX, Y and XY polarizations are for the
                    # interpolation
                    bool_Flags_pre2_X = bool_Flags_X[0:-4]
                    bool_Flags_pre1_X = bool_Flags_X[1:-3]
                    bool_Flags_fol1_X = bool_Flags_X[3:-1]
                    bool_Flags_fol2_X = bool_Flags_X[4:]
                    
                    bool_Flags_pre2_Y = bool_Flags_Y[0:-4]
                    bool_Flags_pre1_Y = bool_Flags_Y[1:-3]
                    bool_Flags_fol1_Y = bool_Flags_Y[3:-1]
                    bool_Flags_fol2_Y = bool_Flags_Y[4:]
                    
                    bool_Flags_pre1_XY = bool_Flags_XY[1:-3]
                    bool_Flags_fol1_XY  =bool_Flags_XY[3:-1]
                    
                    # dentify if the possible polarizations for interpolation are not too far
                    # considering time acquisition
                    bool_t_smos_pre2 = (t_smos_sec_uniq[2:-2]-t_smos_sec_uniq[0:-4]) <= 2.5
                    bool_t_smos_pre1 = (t_smos_sec_uniq[2:-2]-t_smos_sec_uniq[1:-3]) <= 1.3
                    bool_t_smos_fol1 = (t_smos_sec_uniq[3:-1]-t_smos_sec_uniq[2:-2]) <= 1.3
                    bool_t_smos_fol2 = (t_smos_sec_uniq[4:]  -t_smos_sec_uniq[2:-2]) <= 2.5
    
                    # Interpolation possibilities for each polarization X, Y and XY
                    bool_interp_X = np.concatenate(([1<0],[1<0],\
                        np.logical_and(np.logical_or(np.logical_and(bool_Flags_pre2_X,bool_t_smos_pre2), \
                                                     np.logical_and(bool_Flags_pre1_X,bool_t_smos_pre1)),\
                                       np.logical_or(np.logical_and(bool_Flags_fol1_X,bool_t_smos_fol1), \
                                                     np.logical_and(bool_Flags_fol2_X,bool_t_smos_fol2))), \
                                              [1<0],[1<0]),axis=None)                  
                    bool_interp_Y = np.concatenate(([1<0],[1<0],\
                        np.logical_and(np.logical_or(np.logical_and(bool_Flags_pre2_Y,bool_t_smos_pre2), \
                                                     np.logical_and(bool_Flags_pre1_Y,bool_t_smos_pre1)), \
                                       np.logical_or(np.logical_and(bool_Flags_fol1_Y,bool_t_smos_fol1), \
                                                     np.logical_and(bool_Flags_fol2_Y,bool_t_smos_fol2))), \
                                                 [1<0],[1<0]),axis=None)                 
                    bool_interp_XY= np.concatenate(([1<0],[1<0], \
                                      np.logical_and(np.logical_and(bool_Flags_pre1_XY,bool_t_smos_pre1), \
                                                     np.logical_and(bool_Flags_fol1_XY,bool_t_smos_fol1)), \
                                             [1<0],[1<0]),axis=None)
                    
                    if (np.sum(bool_Flags_X) >= 2 and np.sum(bool_Flags_Y) >= 2 and \
                        np.sum(bool_Flags_XY) >= 2):
                        
                        # Interpolation of X, Y and XY where needed
                        I_TBx = np.interp(t_smos_sec_uniq[bool_interp_X], t_smos_sec_uniq[bool_Flags_X],TBxy[bool_Flags_X,0])
                        I_TBy = np.interp(t_smos_sec_uniq[bool_interp_Y],t_smos_sec_uniq[bool_Flags_Y],TBxy[bool_Flags_Y,1])
                        I_TBre = np.interp(t_smos_sec_uniq[bool_interp_XY],t_smos_sec_uniq[bool_Flags_XY],TBxy[bool_Flags_XY,2])
                        I_TBim = np.interp(t_smos_sec_uniq[bool_interp_XY],t_smos_sec_uniq[bool_Flags_XY],TBxy[bool_Flags_XY,3])
                        
                        # Interpolation of the radiometric accuracies
                        I_RATBx = np.interp(t_smos_sec_uniq[bool_interp_X],t_smos_sec_uniq[bool_Flags_X],TBxy[bool_Flags_X,7])
                        I_RATBy = np.interp(t_smos_sec_uniq[bool_interp_Y],t_smos_sec_uniq[bool_Flags_Y],TBxy[bool_Flags_Y,8])
                        I_RATBxy = np.interp(t_smos_sec_uniq[bool_interp_XY],t_smos_sec_uniq[bool_Flags_XY],TBxy[bool_Flags_XY,9])
                        
                        # Saving the interpolated values
                        TBxy[bool_interp_X,0]=I_TBx
                        TBxy[bool_interp_Y,1]=I_TBy
                        TBxy[bool_interp_XY,2]=I_TBre
                        TBxy[bool_interp_XY,3]=I_TBim
                        TBxy[bool_interp_X,7]=I_RATBx
                        TBxy[bool_interp_Y,8]=I_RATBy
                        TBxy[bool_interp_XY,9]=I_RATBxy
                        
                        # ROTATION & ERROR PROPAGATION
                        # Transformation
                        TBhv = np.full([len(Snap_ID_uniq),4],np.nan)
                        RATBhv = np.full([len(Snap_ID_uniq),4],np.nan)
                        
                        # Parallelizing the inversion of MR4 for all angles/Snap_ID - FC 01/11
                        IMR4 = invMR4L(alpha_uniq)
                        TBhv[:,0] = IMR4[:,0,0]*TBxy[:,0] + IMR4[:,0,1]*TBxy[:,1] + IMR4[:,0,2]*TBxy[:,2]*2 - IMR4[:,0,3]*TBxy[:,3]*2
                        TBhv[:,1] = IMR4[:,1,0]*TBxy[:,0] + IMR4[:,1,1]*TBxy[:,1] + IMR4[:,1,2]*TBxy[:,2]*2 - IMR4[:,1,3]*TBxy[:,3]*2
                        TBhv[:,2] = IMR4[:,2,0]*TBxy[:,0] + IMR4[:,2,1]*TBxy[:,1] + IMR4[:,2,2]*TBxy[:,2]*2 - IMR4[:,2,3]*TBxy[:,3]*2
                        TBhv[:,3] = IMR4[:,3,0]*TBxy[:,0] + IMR4[:,3,1]*TBxy[:,1] + IMR4[:,3,2]*TBxy[:,2]*2 - IMR4[:,3,3]*TBxy[:,3]*2
    
                        RATBhv[:,0] = (IMR4[:,0,0]**2*TBxy[:,7]**2 + IMR4[:,0,1]**2*TBxy[:,8]**2 + 4*(IMR4[:,0,2]**2 + IMR4[:,0,3]**2)*TBxy[:,9]**2)**0.5
                        RATBhv[:,1] = (IMR4[:,1,0]**2*TBxy[:,7]**2 + IMR4[:,1,1]**2*TBxy[:,8]**2 + 4*(IMR4[:,1,2]**2 + IMR4[:,1,3]**2)*TBxy[:,9]**2)**0.5
                        RATBhv[:,2] = (IMR4[:,2,0]**2*TBxy[:,7]**2 + IMR4[:,2,1]**2*TBxy[:,8]**2 + 4*(IMR4[:,2,2]**2 + IMR4[:,2,3]**2)*TBxy[:,9]**2)**0.5
                        RATBhv[:,3] = (IMR4[:,3,0]**2*TBxy[:,7]**2 + IMR4[:,3,1]**2*TBxy[:,8]**2 + 4*(IMR4[:,3,2]**2 + IMR4[:,3,3]**2)*TBxy[:,9]**2)**0.5
       
                        # OUTPUT STRUCTURE
                        end_id = start_id + len(Snap_ID_uniq)
                        end_id_list[idx_in_dgg] =  end_id
                        out_lon[start_id:end_id]=np.matlib.repmat(lon[idx_dgg],1,len(Snap_ID_uniq))
                        out_lat[start_id:end_id]=np.matlib.repmat(lat[idx_dgg],1,len(Snap_ID_uniq))
                        out_Tbh[start_id:end_id]=TBhv[:,0]
                        out_Tbv[start_id:end_id]=TBhv[:,1]
                        out_T3[start_id:end_id]=TBhv[:,2]
                        out_T4[start_id:end_id]=TBhv[:,3]
                        out_RAh[start_id:end_id]=RATBhv[:,0]
                        out_RAv[start_id:end_id]=RATBhv[:,1]
                        out_RA3[start_id:end_id]=RATBhv[:,2]
                        out_RA4[start_id:end_id]=RATBhv[:,3]
                        
                        out_Xorg[start_id + np.nonzero(bool_Flags_X == 1)[0]] = 1
                        out_Yorg[start_id + np.nonzero(bool_Flags_Y == 1)[0]] = 1
                        out_RXYorg[start_id + np.nonzero(bool_Flags_RXY == 1)[0]] = 1
                        out_IXYorg[start_id + np.nonzero(bool_Flags_IXY == 1)[0]] = 1
                        
                        out_inc[start_id:end_id] = Theta_uniq
                        out_Az[start_id:end_id] = Az_uniq
                        start_id=end_id 
                        
                    else:
    
                        end_id = start_id + 1
                        end_id_list[idx_in_dgg] = end_id 
                        out_lon[start_id:end_id] = np.nan
                        out_lat[start_id:end_id] = np.nan
                        out_Tbh[start_id:end_id] = np.nan
                        out_Tbv[start_id:end_id] = np.nan
                        out_T3[start_id:end_id] = np.nan
                        out_T4[start_id:end_id] = np.nan
                        out_RAh[start_id:end_id] = np.nan
                        out_RAv[start_id:end_id] = np.nan
                        out_RA3[start_id:end_id] = np.nan
                        out_RA4[start_id:end_id] = np.nan
                        
                        out_Xorg[start_id:end_id] = np.nan
                        out_Yorg[start_id:end_id] =  np.nan
                        out_RXYorg[start_id:end_id] = np.nan
                        out_IXYorg[start_id:end_id] = np.nan
                        out_inc[start_id:end_id] = np.nan
                        out_Az[start_id:end_id] = np.nan
                        start_id = end_id 
                        
            else:
                
                end_id = start_id +1 
                end_id_list[idx_in_dgg] = end_id 
                out_lon[start_id:end_id] = np.nan
                out_lat[start_id:end_id] = np.nan
                out_Tbh[start_id:end_id] = np.nan
                out_Tbv[start_id:end_id] = np.nan
                out_T3[start_id:end_id] = np.nan
                out_T4[start_id:end_id] = np.nan
                out_RAh[start_id:end_id] = np.nan
                out_RAv[start_id:end_id] = np.nan
                out_RA3[start_id:end_id] = np.nan
                out_RA4[start_id:end_id] = np.nan
                
                out_Xorg[start_id:end_id] = np.nan
                out_Yorg[start_id:end_id] =  np.nan
                out_RXYorg[start_id:end_id] = np.nan
                out_IXYorg[start_id:end_id] = np.nan
                out_inc[start_id:end_id] = np.nan
                out_Az[start_id:end_id] = np.nan
                start_id = end_id
                               
        #only output the useful data
        out_lon = out_lon[0:end_id]
        out_lat = out_lat[0:end_id]
        out_Tbh = out_Tbh[0:end_id]
        out_Tbv = out_Tbv[0:end_id]
        out_T3 = out_T3[0:end_id]
        out_T4 = out_T4[0:end_id]
        out_RAh = out_RAh[0:end_id]
        out_RAv = out_RAv[0:end_id]
        out_RA3 = out_RA3[0:end_id]
        out_RA4 = out_RA4[0:end_id]
        
        out_Xorg = out_Xorg[0:end_id]
        out_Yorg = out_Yorg[0:end_id]
        out_RXYorg = out_RXYorg[0:end_id]
        out_IXYorg = out_IXYorg[0:end_id]
        
        out_inc = out_inc[0:end_id]
        out_Az = out_Az[0:end_id]
    
    del dgg_list
    return out_lon, out_lat, out_inc, out_Az, end_id_list, out_Tbh, out_Tbv, \
           out_T3, out_T4, out_RAh, out_RAv, out_RA3, out_RA4, \
           out_Xorg, out_Yorg, out_RXYorg, out_IXYorg
