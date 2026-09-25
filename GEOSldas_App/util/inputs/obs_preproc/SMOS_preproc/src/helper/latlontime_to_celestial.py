import numpy as np
from datetime import datetime
    
def latlontime_to_celestial(lat=None,lon=None,theta=None,phi=None,time=None,*args,**kwargs):

# INPUT:
#  lat   = latitude            [degree] 
#  lon   = longitude           [degree] 
#  theta = incidence angle     [degree]
#  phi   = azimuth             [degree]
#  time  = vector: date + time [Y M D HH MM SS]
    
# OUTPUT:
#  alpha = right ascension     [degree]
#  delta = declination         [degree]
#---------------------------------------------------------------------
    
    rad2deg = 180/np.pi
    theta = theta/rad2deg
    phi = phi/rad2deg
    lat = lat/rad2deg
    lon = lon/rad2deg
    
# find elapsed minutes since Jan 1, 2000 (reference epoch)
#-------------------------------------------

    t1 = datetime(2000,1,1,0,0,0)    
    t2 = time

# time diff between t1 and t2

    UTC = t2-t1                                  
    UTC = UTC.total_seconds()/86400              #in day
    
# convert (lat,lon,az,inc,time) to celestial coordinates
#-------------------------------------------
#(theta, pi+phi) = specular direction points to the location in the celestial sky
#standard astronomical reference local geographic coordinates, 
#elevation el, and azimuth Az
    
    el = np.pi/2 - theta                       #elevation or: 90-theta_in_degrees
    Az = phi                                   #astronom azimuth(0 towards south?)
    
#celestial coordinates = declination (delta), right ascension (alpha)
   
#--- alpha (RA)
# The right ascension is the angle between an object and the location 
# of the vernal equinox (First Point in Aries) measured eastward along 
# the celestial equator in hours, minutes, and seconds of sidereal time. 
# Since the location of the vernal equinox changes due to the precession 
# of the Earth's axis of rotation, coordinates must be given with 
# reference to a date or epoch.
    
    JD2000_DAY = 2451544.5                      #reference epoch, 1 Jan 2000
    JD_CENT = 36525                             #Julian century
    Omega_E = 0.2506846*60.0*24.0               #deg/day; Earth rotation rate; (15 degree/h)
    
    Y = time.year
    M = time.month
    D = time.day
    HH = time.hour
    MM = time.minute
    SS = time.second

    C0 = 1721013.5 + D + (HH + MM/60 + SS/3600)/24
    JD = 367*Y - np.floor(1.75*(Y + np.floor((M + 9)/12))) + np.floor(275*M/9) + C0  #julian date
    
    if np.abs(JD - (UTC + JD2000_DAY)) > 1:
        print('is the JD calculated correctly?')
    
    U0 = (JD - JD2000_DAY) / JD_CENT
    Th_G0 = 100.46062 + 36000.77*U0 + 0.000388*U0**2 - 2.6*10**(-8)*U0**3    # degree
    Th_G0 = np.mod(Th_G0,360)

    Th_L = Th_G0 + Omega_E*np.mod(UTC,1) + lon*rad2deg   # degree; local sideral time
    Th_L = np.mod(Th_L,360)

    tmp_H = np.sin(Az) / (np.tan(el)*np.cos(lat) + np.cos(Az)*np.sin(lat))
    H = np.arctan(tmp_H)*rad2deg

    alpha = np.mod(Th_L-H, 360)
    
#--- delta (DEC)
# The declination of an object is its angle in degrees, minutes, and seconds of arc
# above or below the celestial equator.
    
    tmp_delta = np.sin(lat)*np.sin(el) - np.cos(lat)*np.cos(el)*np.cos(Az)
    delta = np.arcsin(tmp_delta) * rad2deg
    
    return alpha, delta, Th_G0, UTC


