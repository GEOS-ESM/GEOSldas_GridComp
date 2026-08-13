import os
import numpy as np
import xml.etree.ElementTree as ET
from src.readwrite import read_aux_gal_sm
from src.helper import latlontime_to_celestial


#-------------------------------------------------------------------------
# Remove galactic and atmospheric contributions.
# Correct SMOS Tb in order to mimic SMAP.
#-------------------------------------------------------------------------
    
# INPUT: Tb_H, Tb_V  Tb                           [K] 
#        lat, lon    latitude, longitude          [degree]
#        angle, phi  incidence, azimuth angle     [degree]
#        time        time                         [Y M D HH MM SS]
#        T_air       surface air temperature      [oC]
#        T_surf      surface temperature          [oC]
#        P_surf      surface pressure             [mbar] 
#        V_surf      surface water vapour density [g/m^3]
#        C           efficiency                   [-]
    
#=========================================================================
def sub2ind(array_shape, rows, cols):
    ind = rows*array_shape[1] + cols
    ind[ind < 0] = -1
    ind[ind >= array_shape[0]*array_shape[1]] = -1
    return ind

def ind2sub(array_shape, ind):
    ind[ind < 0] = -1
    ind[ind >= array_shape[0]*array_shape[1]] = -1
    rows = (ind.astype('int') / array_shape[1])
    cols = ind % array_shape[1]
    return (rows, cols)

def gal_atm_correction(Tb_H=None,Tb_V=None,lat=None,lon=None,angle=None, \
                           phi=None,time=None,T_air=None,T_surf=None, \
                        P_surf=None,V_surf=None,C=None,*args,**kwargs):

    d2r=np.pi / 180
    C_2_K=273.15

# 0) reduce to good data:
#----------------------------
    
    length_org=Tb_H.size
    limited_ind = np.nonzero( \
                   np.logical_not(np.isnan(lat)) &  \
                   np.logical_not(np.isnan(lon)) &   \
                   np.logical_not(np.isnan(angle)) & \
                   np.logical_not(np.isnan(phi)) & \
                   np.logical_not(np.isnan(T_surf)) & \
                   np.logical_not(np.isnan(P_surf)) & \
                   np.logical_not(np.isnan(V_surf)))[0]    

    print('galactic and atmospheric correction on ', str(len(limited_ind)),' of ',length_org)
    Tb_ap_H_out = np.full([length_org],np.nan)
    Tb_ap_V_out = np.full([length_org],np.nan)
    Tb_BOA_H_out = np.full([length_org],np.nan)
    Tb_BOA_V_out = np.full([length_org],np.nan)

    Tb_H=Tb_H[limited_ind]
    Tb_V=Tb_V[limited_ind]
    lat=lat[limited_ind]
    lon=lon[limited_ind]
    angle=angle[limited_ind]
    phi=phi[limited_ind]
    T_air=T_air[limited_ind]
    T_surf=T_surf[limited_ind]
    P_surf=P_surf[limited_ind]
    V_surf=V_surf[limited_ind]

    if Tb_H.size > 0:
        # 1) Remove galactic terms:

        # Get Tc+Tgal from SMOS AUX data

        # v724
        current_dir = os.path.dirname(os.path.abspath(__file__))
        data_file=os.path.join(current_dir,'..','..','data','SM_OPER_AUX_GAL_SM_20050101T000000_20500101T000000_001_003_3')

        print('reading galaxy info')
        p=data_file+'/'+data_file[-60:]+'.HDR'
        doc = ET.parse(p)
        root = doc.getroot()
        
        for tmp in root.iter("{http://213.170.46.150/smos/schemas}Min_RA"):
            Min_ra = float(tmp.text)
            del tmp
 
        for tmp in root.iter("{http://213.170.46.150/smos/schemas}Max_RA"):
            Max_ra = float(tmp.text)    
            del tmp
        
        for tmp in root.iter("{http://213.170.46.150/smos/schemas}Min_DEC"):
            Min_dec = float(tmp.text)    
            del tmp
            
        for tmp in root.iter("{http://213.170.46.150/smos/schemas}Max_DEC"):
            Max_dec = float(tmp.text)    
            del tmp
            
        for tmp in root.iter("{http://213.170.46.150/smos/schemas}DELTA_RA"):
            D_ra = float(tmp.text)    
            del tmp

        for tmp in root.iter("{http://213.170.46.150/smos/schemas}DELTA_DEC"):
            D_dec = float(tmp.text)    
            del tmp

        TB_Sky_H,TB_Sky_V=read_aux_gal_sm(data_file+'/'+data_file[-60:]+'.DBL')

#First row corresponds to largest declination.
#We need to reverse the row order from the smallest to the largest
#declination. This was *not* done in the orginal code used for the GRSL
#paper [bug in paper].

        TB_Sky_H=TB_Sky_H[::-1,:]
        TB_Sky_V=TB_Sky_V[::-1,:]

        alpha,delta,Th_G0,UTC=latlontime_to_celestial(lat,lon,angle,phi,time)

        idx = np.ravel_multi_index(  \
            [np.round((delta-Min_dec)/D_dec).astype('int'), \
             np.round((alpha-Min_ra)/D_ra).astype('int')], \
             TB_Sky_H.shape,  mode='raise', order='F')
            
        TcTgal_H=TB_Sky_H.flatten('F')[idx]
        TcTgal_V=TB_Sky_V.flatten('F')[idx]

        #print('min-max Tsky (H) = ' + str(TcTgal_H.min()) + ' ' + str(TcTgal_H.max()))
        #print('min-max Tsky (V) = ' + str(TcTgal_V.min()) + ' ' + str(TcTgal_V.max()))
        
        # Estimate epsilon
        eps_H = C * Tb_H / (T_surf + C_2_K)
        eps_V = C * Tb_V / (T_surf + C_2_K)
        
        #print('reset eps (H) = ' + str(np.nansum(eps_H > 1)) + ' out of ' + str(eps_H.size))
        #print('reset eps (V) = ' + str(np.nansum(eps_V > 1)) + ' out of ' + str(eps_V.size))
        
        eps_H[eps_H > 1]=1
        eps_V[eps_V > 1]=1

        #print('min-max eps (H) = ' + str(eps_H.min()) + ' ' + str(eps_H.max()))
        #print('min-max eps (V) = ' + str(eps_V.min()) + ' ' + str(eps_V.max()))
        
        # Estimate atmospheric loss factor:
        # coefficients according to L1B ATBD (slightly different from Peng et al., 2013)
        L=1.00938 - 2.9626*10**(-5)*T_air + 1.6521*10**(-5)*(P_surf-900) + 1.0712*10**(-5)*V_surf

        #print('min-max T_air = ' + str(T_air.min()) + ' ' + str(T_air.max()))
        #print('min-max T_surf = ' + str(T_surf.min())+ ' ' + str(T_surf.max()))
        #print('min-max V_surf = ' + str(V_surf.min())+ ' ' + str(V_surf.max()))

        ang_corr=np.cos(40.0*d2r) * 1/np.cos(angle*d2r)
        L = L**ang_corr

        #print('min-max angle = ' +str(angle.min()) + ' ' + str(angle.max()))
        #print('min-max atm loss L = ' + str(L.min()) + ' ' +str(L.max()))
        
        # Remove galactic and cosmic radiation
        H_gal_corr = TcTgal_H * (1 - eps_H) / (L*L)
        V_gal_corr = TcTgal_V * (1 - eps_V) / (L*L)

        Tb_ap_H = Tb_H - H_gal_corr
        Tb_ap_V = Tb_V - V_gal_corr

        #print('min-max gal corr (H) = -(' + str(H_gal_corr.min()) + ' ' + str(H_gal_corr.max())+')')
        #print('min-max gal corr (V) = -(' + str(V_gal_corr.min()) + ' ' + str(V_gal_corr.max())+')')
        
        # 2) Remove atmospheric effects:

        # Estimate upwelling Tb:
        # T_up, temporary, coefficients according to L1B ATBD (slightly different from Peng et al., 2013)
        T_up=2.3058 - 3.2735*10**(-3)*T_air + 4.233*10**(-3)*(P_surf - 900) + 1.4472*10**(-3)*V_surf

        # angular function to correct T_up, coefficients according to Peng et al., 2013
        f_ang=np.nan*angle

        ind = np.nonzero(angle < 20)[0]

        if ind.size:
            f_ang[ind] = 1.2855*10**(-4)*angle[ind]**2 - 1.3361*10**(-4)*angle[ind] + 0.7625

        ind=np.nonzero(np.logical_and( angle >= 20, angle <= 60))[0]

        if ind.size:
            f_ang[ind] = 8.2724*10**(-6)*angle[ind]**3 - 5.7129*10**(-4)*angle[ind]**2 + 2.0411*10**(-2)*angle[ind] + 0.5655

        ind=np.nonzero(np.logical_and(angle > 60,angle <= 70))[0]

        if ind.size:
            f_ang[ind]=2.4189*10**(-3)*angle[ind]**2 - 0.2458*angle[ind]+ 7.5624

        # get final T_up
        T_up = T_up*f_ang

        #print('min-max T_up = ' + str(T_up.min()) + ' ' + str(T_up.max()))
        
        # Remove atmospheric correction, using Eq. below Eq. 5.55 in L1B ATBD
        Tb_BOA_H = (T_surf + C_2_K)*(Tb_ap_H*L - (1 + L)*T_up) / (T_surf + C_2_K - T_up)
        Tb_BOA_V = (T_surf + C_2_K)*(Tb_ap_V*L - (1 + L)*T_up) / (T_surf + C_2_K - T_up)

        H_atm_corr=Tb_BOA_H - Tb_ap_H
        V_atm_corr=Tb_BOA_V - Tb_ap_V

        # prevent an increase in Tb 
        # an increase would happen occasionally, 
        # especically for V-pol and higher incidence angles,
        # when T_up becomes larger and Tb_BOA_V > Tb_ap_V
        #print('reset atm corr (H) = ' + str(np.nansum(H_atm_corr > 0)) + ' out of ' + str(H_atm_corr.size))
        #print('reset atm corr (V) = ' + str(np.nansum(V_atm_corr > 0)) + ' out of ' + str(V_atm_corr.size))
        
        H_atm_corr[H_atm_corr > 0]=0
        V_atm_corr[V_atm_corr > 0]=0

        Tb_BOA_H=Tb_ap_H + H_atm_corr
        Tb_BOA_V=Tb_ap_V + V_atm_corr

        #print('min-max atm corr (H) = ' + str(H_atm_corr.min()) + ' ' + str(H_atm_corr.max()))
        #print('min-max atm corr (V) = ' + str(V_atm_corr.min()) + ' ' + str(V_atm_corr.max()))
        # 3) Expand to original length:
#--------------------------
        Tb_ap_H_out[limited_ind]=Tb_ap_H
        Tb_ap_V_out[limited_ind]=Tb_ap_V

        Tb_BOA_H_out[limited_ind]=Tb_BOA_H
        Tb_BOA_V_out[limited_ind]=Tb_BOA_V
        
    return Tb_ap_H_out, Tb_ap_V_out, Tb_BOA_H_out, Tb_BOA_V_out
    
