import numpy as np
    
#  SMAPEASE2INVERSE  The principal function is to perform inverse transformation
#                    from (row,col)'s to (lat,lon)'s for a set of nested EASE
#                    grids defined at 1, 3, 9, and 36km grid resolutions.  These
#                    grids are all based on the EASE-Grid 2.0 specification (WGS84
#                    ellipsoid).
    
#  SYNTAX            [lat,lon] = smapease2inverse(row,col,gridid)
    
#                    where gridid is a 3-character string enclosed in single
#                    quotes, in the form of {M|N|S}{01,03,09,36}.  This subroutine
#                    accepts vector inputs and produce vector outputs.
    
#  HISTORY           This subroutine was adapted from the offical EASE-Grid-2.0
#                    conversion utilities (written in IDL) developed by the
#                    NSIDC.
    
#                    Note that in NSIDC's original implementation, (row,col) are
#                    zero-based.  In other words, the first cell is (0,0) and the
#                    last cell is (N-1,M-1), where N and M are the row and column
#                    dimensions of the array.  In this MATLAB implementation, the
#                    same convention is used.  In other words, the end point of
#                    the first cell is located at (r,c) = (-0.5,-0.5) whereas the
#                    end point of the last cell is located at (r,c) = (14615.5,
#                    34703.5).  Thus,
    
#                   [lat,lon] = smapease2inverse(-0.5,-0.5,'M01') returns:
#                    lat = 85.044566407398861
#                    lon = 1.799999999999994e+02
    
#                   [lat,lon] = smapease2inverse(14615.5,34703.5,'M01') returns:
#                    lat = -85.044566407398861
#                    lon = -1.799999999999994e+02
    
#                    The polar grids, on the other hand, are more complete in
#                    terms of latitude coverage:
    
#                   [lat,lon] = smapease2inverse(8999,8999,'N01')
#                    lat = 89.993669248945238 
#                    lon = -135
#                   [lat,lon] = smapease2inverse(9000,9000,'N01')
#                    lat = 89.993669248945238
#                    lon = 45
    
#                   [lat,lon] = smapease2inverse(8999,8999,'S01')
#                    lat = -89.993669248945238 
#                    lon = -45
#                   [lat,lon] = smapease2inverse(9000,9000,'S01')
#                    lat = -89.993669248945238
#                    lon = 135
    
#  UPDATE            North/south polar projections were added. (03/2012)
    
#  REFERENCE         Brodzik, M. J., B. Billingsley, T. Haran, B. Raup, and M. H.
#                    Savoie (2012): EASE-Grid 2.0: Incremental but Significant
#                    Improvements for Earth-Gridded Data Sets. ISPRS International
#                    Journal of Geo-Information, vol. 1, no. 1, pp. 32-45,
#                    http://www.mdpi.com/2220-9964/1/1/32/
#  
#  Steven Chan, 11/2011
#  Email: steven.k.chan@jpl.nasa.gov
    
def EASEv2_ind2latlon(row=None,col=None,gridid=None,*args,**kwargs):

    # Constants returned by EASE2_GRID_INFO.PRO
    projection=gridid[0]
    if 'M36' == gridid:
        map_scale_m=36032.220840584
        cols=964
        rows=406
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'M09' == gridid:
        map_scale_m=9008.055210146
        cols=3856
        rows=1624
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'M03' == gridid:
        map_scale_m=3002.6850700487
        cols=11568
        rows=4872
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'M01' == gridid:
        map_scale_m=1000.89502334956
        cols=34704
        rows=14616
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'N36' == gridid:
        map_scale_m=36000.0
        cols=500
        rows=500
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'N09' == gridid:
        map_scale_m=9000.0
        cols=2000
        rows=2000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'N03' == gridid:
        map_scale_m=3000.0
        cols=6000
        rows=6000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'N01' == gridid:
        map_scale_m=1000.0
        cols=18000
        rows=18000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'S36' == gridid:
        map_scale_m=36000.0
        cols=500
        rows=500
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'S09' == gridid:
        map_scale_m=9000.0
        cols=2000
        rows=2000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'S03' == gridid:
        map_scale_m=3000.0
        cols=6000
        rows=6000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'S01' == gridid:
        map_scale_m=1000.0
        cols=18000
        rows=18000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    else:
        print('ERROR: Incompatible grid specification.')
    
    # Constants returned by EASE2_MAP_INFO.PRO
    map_equatorial_radius_m=6378137.0
    map_eccentricity=0.081819190843
    e2=map_eccentricity ** 2
    if 'M' == projection:
        map_reference_latitude=0.0
        map_reference_longitude=0.0
        map_second_reference_latitude=30.0
        sin_phi1=np.sin(map_second_reference_latitude*np.pi / 180)
        cos_phi1=np.cos(map_second_reference_latitude*np.pi / 180)
        kz=cos_phi1 / np.sqrt(1.0 - e2*sin_phi1*sin_phi1)
    elif 'N' == projection:
        map_reference_latitude=90.0
        map_reference_longitude=0.0
    elif 'S' == projection:
        map_reference_latitude=- 90.0
        map_reference_longitude=0.0
    
    # Selected calculations inside WGS84_INVERSE.PRO
    x= (col - r0) * map_scale_m
    y= (s0 - row) * map_scale_m
    # Selected calculations inside WGS84_INVERSE_XY.PRO
    e4=map_eccentricity ** 4
    e6=map_eccentricity ** 6
    qp= (1.0 - e2)*((1.0 / (1.0 - e2)) - (1.0 / (2.0*map_eccentricity))*np.log((1.0 - map_eccentricity) / (1.0 + map_eccentricity)))
    if 'M' == projection:
        beta= np.arcsin(2.0 * y * kz / (map_equatorial_radius_m * qp))
        lam= x / (map_equatorial_radius_m*kz)
    elif 'N' == projection:
        rho=np.sqrt(x ** 2 + y ** 2)
        beta=np.arcsin(1.0 - (rho ** 2) / (qp * map_equatorial_radius_m ** 2))
        lam=np.arctan2(x,- y)
    elif 'S' == projection:
        rho=np.sqrt(x ** 2 + y ** 2)
        beta=- np.arcsin(1.0 - (rho ** 2) / (qp * map_equatorial_radius_m ** 2))
        lam=np.arctan2(x,y)
    
    phi= beta + (((e2 / 3.0) + ((31.0 / 180.0) * e4) + ((517.0 / 5040.0) * e6)) * np.sin(2.0 * beta)) + \
        ((((23.0 / 360.0) * e4) + ((251.0 / 3780.0) * e6)) * np.sin(4.0*beta)) + (((761.0 / 45360.0) * e6) * np.sin(6.0 * beta))
    lat= phi * 180.0 / np.pi
    lon= map_reference_longitude + (lam * 180.0 / np.pi)

    lat  = np.atleast_1d(np.asarray(lat,  dtype=np.float64))
    lon  = np.atleast_1d(np.asarray(lon,  dtype=np.float64))
    
    msk1= np.where(lon < - 180.0)
    lon[msk1]= lon[msk1] + 360.0
    msk2= np.where(lon >   180.0)
    lon[msk2]=lon[msk2] - 360.0
    if 'N' == projection:
        idx=np.where(lat < 0.0)
        lat[idx]= float('nan')
        lon[idx]= float('nan')
    elif 'S' == projection:
        idx= np.where(lat > 0.0)
        lat[idx]= float('nan')
        lon[idx]= float('nan')
            
    return lat, lon        

# --------------------------------------------------------------------------------------------------------

#  SMAPEASE2FORWARD  The principal function is to perform forward transformation
#                    from (lat,lon)'s to (row,col)'s for a set of nested EASE
#                    grids defined at 1, 3, 9, and 36km grid resolutions.  These
#                    grids are all based on the EASE-Grid 2.0 specification (WGS84
#                    ellipsoid).
    
#  SYNTAX            [row,col] = smapease2forward(lat,lon,gridid)
    
#                    where gridid is a 3-character string enclosed in single
#                    quotes, in the form of {M|N|S}{01,03,09,36}.  This subroutine
#                    accepts vector inputs and produce vector outputs.
    
#  HISTORY           This subroutine was adapted from the offical EASE-Grid-2.0
#                    conversion utilities (written in IDL) developed by the
#                    NSIDC.
    
#                    Note that in NSIDC's original implementation, (row,col) are
#                    zero-based.  In other words, the first cell is (0,0) and the
#                    last cell is (N-1,M-1), where N and M are the row and column
#                    dimensions of the array.  In this MATLAB implementation, the
#                    same convention is used.  In other words, the end point of
#                    the first cell is located at (r,c) = (-0.5,-0.5) whereas the
#                    end point of the last cell is located at (r,c) = (14615.5,
#                    34703.5).  Thus,
    
#                   [lat,lon] = smapease2inverse(-0.5,-0.5,'M01') returns:
#                    lat = 85.044566407398861
#                    lon = 1.799999999999994e+02
    
#                   [lat,lon] = smapease2inverse(14615.5,34703.5,'M01') returns:
#                    lat = -85.044566407398861
#                    lon = -1.799999999999994e+02
    
#                    The polar grids, on the other hand, are more complete in
#                    terms of latitude coverage:
    
#                   [lat,lon] = smapease2inverse(8999,8999,'N01')
#                    lat = 89.993669248945238 
#                    lon = -135
#                   [lat,lon] = smapease2inverse(9000,9000,'N01')
#                    lat = 89.993669248945238
#                    lon = 45
    
#                   [lat,lon] = smapease2inverse(8999,8999,'S01')
#                    lat = -89.993669248945238 
#                    lon = -45
#                   [lat,lon] = smapease2inverse(9000,9000,'S01')
#                    lat = -89.993669248945238
#                    lon = 135
    
#  UPDATE            North/south polar projections were added. (03/2012)
    
#  REFERENCE         Brodzik, M. J., B. Billingsley, T. Haran, B. Raup, and M. H.
#                    Savoie (2012): EASE-Grid 2.0: Incremental but Significant
#                    Improvements for Earth-Gridded Data Sets. ISPRS International
#                    Journal of Geo-Information, vol. 1, no. 1, pp. 32-45,
#                    http://www.mdpi.com/2220-9964/1/1/32/
#  
#  Steven Chan, 11/2011
#  Email: steven.k.chan@jpl.nasa.gov
    
def EASEv2_latlon2ind(lat=None,lon=None,gridid=None,*args,**kwargs):
    lat  = np.atleast_1d(np.asarray(lat,  dtype=np.float64))
    lon  = np.atleast_1d(np.asarray(lon,  dtype=np.float64))
    # Constants returned by EASE2_GRID_INFO.PRO
    projection=gridid[0]
    if 'M36' == gridid:
        map_scale_m=36032.220840584
        cols=964
        rows=406
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'M09' == gridid:
        map_scale_m=9008.055210146
        cols=3856
        rows=1624
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'M03' == gridid:
        map_scale_m=3002.6850700487
        cols=11568
        rows=4872
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'M01' == gridid:
        map_scale_m=1000.89502334956
        cols=34704
        rows=14616
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'N36' == gridid:
        map_scale_m=36000.0
        cols=500
        rows=500
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'N09' == gridid:
        map_scale_m=9000.0
        cols=2000
        rows=2000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'N03' == gridid:
        map_scale_m=3000.0
        cols=6000
        rows=6000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'N01' == gridid:
        map_scale_m=1000.0
        cols=18000
        rows=18000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'S36' == gridid:
        map_scale_m=36000.0
        cols=500
        rows=500
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'S09' == gridid:
        map_scale_m=9000.0
        cols=2000
        rows=2000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'S03' == gridid:
        map_scale_m=3000.0
        cols=6000
        rows=6000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    elif 'S01' == gridid:
        map_scale_m=1000.0
        cols=18000
        rows=18000
        r0=(cols - 1) / 2
        s0=(rows - 1) / 2
    else:
        print('ERROR: Incompatible grid specification.')
    
    # Constants returned by EASE2_MAP_INFO.PRO
    epsilon=1e-06
    map_equatorial_radius_m=6378137.0
    map_eccentricity=0.081819190843
    e2=map_eccentricity ** 2
    if 'M' == projection:
        map_reference_latitude=0.0
        map_reference_longitude=0.0
        map_second_reference_latitude=30.0
        sin_phi1=np.sin(map_second_reference_latitude*np.pi / 180)
        cos_phi1=np.cos(map_second_reference_latitude*np.pi / 180)
        kz=cos_phi1 / np.sqrt(1.0 - e2*sin_phi1*sin_phi1)
    elif 'N' == projection:
        map_reference_latitude=90.0
        map_reference_longitude=0.0
    elif 'S' == projection:
        map_reference_latitude=-90.0
        map_reference_longitude=0.0
        
    # Selected calculations inside WGS84_CONVERT.PRO and WGS84_CONVERT_XY.PRO
    dlon=lon - map_reference_longitude
    msk1= dlon < -180.0
    dlon[msk1]=dlon[msk1] + 360.0
    msk2= dlon >  180.0
    dlon[msk2]=dlon[msk2] - 360.0
    phi= lat*np.pi / 180.0
    lam= dlon*np.pi / 180.0
    sin_phi=np.sin(phi)
    q=(1.0 - e2)*((sin_phi / (1.0 - e2*sin_phi*sin_phi)) - (1.0 / (2.0*map_eccentricity))*np.log((1.0 - map_eccentricity * sin_phi) / (1.0 + map_eccentricity * sin_phi)))
    qp=1.0 - ((1.0 - e2) / (2.0*map_eccentricity)*np.log((1.0 - map_eccentricity) / (1.0 + map_eccentricity)))
    if 'M' == projection:
        x= map_equatorial_radius_m * kz * lam
        y= (map_equatorial_radius_m * q) / (2.0*kz)
    elif 'N' == projection:
        tmp=qp - q
        tmp[np.absolute(tmp) < epsilon]=0.0
        rho= map_equatorial_radius_m*np.sqrt(tmp)
        x= rho * np.sin(lam)
        y= -rho * np.cos(lam)
    elif 'S' == projection:
        tmp=qp + q
        tmp[np.absolute(tmp) < epsilon]=0.0
        rho= map_equatorial_radius_m * np.sqrt(tmp)
        x= rho * np.sin(lam)
        y= rho * np.cos(lam)

    row=s0 - (y / map_scale_m)
    col=r0 + (x / map_scale_m)

    if 'N' == projection:
        idx=(lat < 0.0)
        row[idx]=np.nan
        col[idx]=np.nan
    elif 'S' == projection:
        idx=(lat > 0.0)
        row[idx]=np.nan
        col[idx]=np.nan
    
    # restrict indices to (zero-based) min and max indices
    # (needed because of round-off error - GDL, 07 Jun 2013)
    
    row=np.maximum(0,np.minimum(rows - 1,np.round(row)))
    col=np.maximum(0,np.minimum(cols - 1,np.round(col)))

    # assign max index if nan (unclear why. this is based on the matlab code)
    row[(np.isnan(row))] = rows - 1
    col[(np.isnan(col))] = cols - 1
    return row.astype('int'),col.astype('int')
    

# ============= EOF ========================================================================================
