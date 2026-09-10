from netCDF4 import Dataset 
import numpy as np

def isKthBitSet(number, kth): 
    if number & (1 << (kth - 1)): 
        return True 
    else: 
        return False 

def bitget(number, kth):
    return (number >> kth-1) & 1

def ismember(a, b):
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    # None can be replaced by any other "not in b" value
    return [bind.get(itm, None) for itm in a]  

def read_SCLF1C_nc(fname=None,*args,**kwargs):
    
    fin = Dataset(fname,'r')
    
    lat = fin.variables['Grid_Point_Latitude'][:]
    lon = fin.variables['Grid_Point_Longitude'][:]
    
    mask_grid = fin.variables['Grid_Point_Mask'][:]
    
    mask_ok = np.logical_and(bitget(mask_grid,2), np.logical_not(bitget(mask_grid,6)))
     
    Grid_ID = fin.variables['Grid_Point_ID'][:]
    flag_grid = fin.variables['Flags'][:]
    
    RA = fin.variables['Radiometric_Accuracy_of_Pixel'][:]

    Theta = 1.*fin.variables['Incidence_Angle'][:]
    Az = fin.variables['Azimuth_Angle'][:]
    Fa = fin.variables['Faraday_Rotation_Angle'][:]
    Ge = fin.variables['Geometric_Rotation_Angle'][:]
    
    Seconds = fin.variables['Seconds'][:]
    Microseconds = fin.variables['Microseconds'][:]

    BT_Value_Real = fin.variables['BT_Value_Real'][:] 
    BT_Value_Imag = fin.variables['BT_Value_Imag'][:]
   
    # v620 
    #mask_nonRFI = np.all((np.logical_not(bitget(flag_grid,12)),np.logical_not(bitget(flag_grid,16)), \
    #                             bitget(flag_grid,11), BT_Value_Real < 325., BT_Value_Imag < 325.), axis=0)
    # v724 
    mask_nonRFI = np.all((np.logical_not(bitget(flag_grid,12)),np.logical_not(bitget(flag_grid,7)), \
                                 bitget(flag_grid,11), BT_Value_Real < 325., BT_Value_Imag < 325.), axis=0)
    
    TB_real = BT_Value_Real
    TB_imag = BT_Value_Imag
    
    TB_real[np.logical_not(mask_nonRFI)] = np.nan
    TB_imag[np.logical_not(mask_nonRFI)] = np.nan
      
    RA[np.logical_not(mask_nonRFI)] = np.nan
    Theta[np.logical_not(mask_nonRFI)] = np.nan
    Az[np.logical_not(mask_nonRFI)] = np.nan
    Fa[np.logical_not(mask_nonRFI)] = np.nan
    Ge[np.logical_not(mask_nonRFI)] = np.nan

    Snap_ID = fin.variables['Snapshot_ID'][:]
    Snap_ID_of_Pixel = fin.variables['Snapshot_ID_of_Pixel'][:]
    
    idx_nonRFI = np.nonzero(mask_nonRFI)
    t_smos_sec = np.full(TB_real.shape,np.nan)
    loc_snap = ismember(Snap_ID_of_Pixel[idx_nonRFI], Snap_ID)
    t_smos_sec[idx_nonRFI] = Seconds[loc_snap] + Microseconds[loc_snap] / 1.e6
    
    flag_16 = bitget(flag_grid,1)
    flag_15 = bitget(flag_grid,2)
    
    BT_count = np.sum(np.logical_not(np.isnan(TB_real)).astype(int),axis=1)
        
    return TB_real, TB_imag, RA, Theta, Az, Fa, Ge, Snap_ID_of_Pixel, t_smos_sec, \
            lat, lon, flag_15, flag_16, mask_ok, Grid_ID, BT_count
            
           
