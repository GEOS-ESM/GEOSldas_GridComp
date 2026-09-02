import numpy as np 
from netCDF4 import Dataset
from EASEv2 import EASEv2_latlon2ind

def read_GEOSIT(date_time=None,vegcls_grid=None,EASEv2_M36_to_M2_grid=None,lat=None,lon=None,GEOSIT_path=None):
    
    # Read fields from a GEOSIT 

    # INPUT: 
    # date_time  =                #structure with year, month, day, hour,...
    # vegcls_grid           =     #2-d vegcls on EASEv2_M36 output grid
    # EASEv2_M36_to_M2_grid =     #2-d EASEv2_M36 grid where each index 
    #                              refers to the nearest GEOSIT LL grid cell
    # lat     = [ : ];            #EASEv2 M36 indices for each element in swath
    # lon     = [ : ];
        
    # OUTPUT:
    # Tair, Ts, Vs, Ps
    # - in units needed for atm and gal correction
    # - for locations within SMOS swath only
    #-----------------------------------------------------------------
    
    C_2_K=273.15
    
    if (date_time.year < 2008):
        print('did not look into the stream sequence for older years than 2010')
    else:
        if (date_time.year < 2018):
            stream='d5294_geosit_jan08'
        else:
            stream='d5294_geosit_jan18'                
 
    #GEOS-IT:
    
    path=GEOSIT_path+stream+'/diag/Y'+date_time.strftime("%Y")+'/M'+date_time.strftime("%m")+'/'

    date_string=date_time.strftime("%Y-%m-%dT%H")

    filename= Dataset(path+'/'+stream+'.lnd_tavg_1hr_glo_L576x361_slv.'+date_string+'30Z.nc4','r')
    
    Tp = filename.variables['TSURF']

    Ts = filename.variables['TSOIL1']

    filename = Dataset(path+'/'+stream+'.slv_tavg_1hr_glo_L576x361_slv.'+date_string+'30Z.nc4','r')

    Ps = filename.variables['PS']

    SH = filename.variables['QV2M']

    Tair = filename.variables['T2M']
    
    #print('Read GEOSIT files - '+stream+'.lnd_tavg_1hr_glo_L576x361_slv.'+date_string+'30Z.nc4')
    # Interpolate, limit to swath
    #-----------------------------------------------------------------
    
    row_ind,col_ind = EASEv2_latlon2ind(lat,lon,'M36')
    idx_M36 = np.ravel_multi_index([row_ind,col_ind],EASEv2_M36_to_M2_grid.shape,  mode='raise', order='F')
    #idx = sub2ind(EASEv2_M36_to_M2_grid.shape,row_ind,col_ind)

    idx = EASEv2_M36_to_M2_grid.ravel(order='F')[idx_M36]
    idx = idx -1 # change to 0-based index
 
    #--- 0) Surface air temperature in [^o C]-------------------------
    
    Tair = Tair[0,:,:]- C_2_K
    Tair = Tair.flatten('F')[idx]
    #--- 1) Surface temperature in [^o C]-----------------------------
    
    Ts = Ts[0,:,:] - C_2_K    
    Tp = Tp[0,:,:] - C_2_K
    
    Ts=Ts.flatten('F')[idx]
    Tp=Tp.flatten('F')[idx]
    Ts[vegcls_grid.flatten('F')[idx_M36] != 1] = \
    (Tp[vegcls_grid.flatten('F')[idx_M36] != 1] + \
     Ts[vegcls_grid.flatten('F')[idx_M36] != 1]) / 2.0
        
    if any(Ts > 70) or any(Ts < - 80):
        print('WARNING: Surface temperature out of range')
    
    #print('Surface temperature min '+str(Ts.min())+' max '+str(Ts.max())+' [^o C]')
    #--- 2) Surface water vapour density [g/m3]------------------------
    
    # SH = specific humidity, ratio of the water vapor content of the 
    # mixture to the total air content on a mass basis.
    # Assuming water density at 1 kg/m^3:
    
    Vs = SH[0,:,:] * 10**3
    
    Vs=Vs.flatten('F')[idx]
    
    if any(Vs > 80):
        print('WARNING: Surface vapour pressure density out of range')
    
    if any(Vs < 0):
        print('WARNING: min vapour pressure density will be reset. min: '+min(Vs))
        Vs[Vs < 0] = 0
    
    #print('Surface water vapour density min '+str(Vs.min())+' max '+str(Vs.max())+' [g/m3]')
    #--- 3) Surface pressure [mbar]------------------------------------
    
    Ps = Ps[0,:,:] * 10**(-5) * 10**3    #  [Pa] --> [mbar]

    
    Ps = Ps.flatten('F')[idx]

    if any(Ps > 1200) or any(Ps < 400):
        print('WARNING: Surface pressure out of range')
    
    #print('Surface pressure min '+str(Ps.min())+' max '+str(Ps.max())+' [mbar]')
    #-----------------------------------------------------------------
    return Tair, Ts, Vs, Ps
    
