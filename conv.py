
# A collection of subroutines to do unit conversions (and a few other things).

import math, time, string
from cons import *

# abetuv and defuvp work together to give a distance on the earth between two lat/lon
#  points.  To use call defuvp twice (once with each set of lat/lon values).  Then
#  call abetuv with the Unit_vector returned from defuvp in A and B for each set of
#  lat/lon coordinates.  It returns the distance between them in meters.

def abetuv(A,B):
#    pi = 3.141592653589793
    cos_pi_over_4 = 0.707106781

    # Use dot product for large angles
    dot_prod = A[0]*B[0]+A[1]*B[1]+A[2]*B[2]
    if (abs(dot_prod) < cos_pi_over_4):
        #Theta[0] = math.acos(dot_prod)
        theta = math.acos(dot_prod)
        
    # Use cross product for small angles
    else:
        cross_prod1 = A[1]*B[2] - A[2]*B[1]
        cross_prod2 = A[2]*B[0] - A[0]*B[2]
        cross_prod3 = A[0]*B[1] - A[1]*B[0]
        cmag = math.sqrt(math.pow(cross_prod1,2) + math.pow(cross_prod2,2) + math.pow(cross_prod3,2))
        #Theta[0] = math.asin(cmag)
        theta = math.asin(cmag)
        if (dot_prod < 0.0):
            #Theta[0] = pi - Theta[0]
            theta = pi - theta
    return (re * theta)
        
def defuvp(lat, lon, Unit_vector):
#    degrees_per_radian = 57.29577951308232
    degrees_per_radian = 180.0/pi

    latr = lat / degrees_per_radian
    lonr = lon / degrees_per_radian
    cos_lat = math.cos(latr)
    Unit_vector[0] = cos_lat * math.cos(lonr)
    Unit_vector[1] = cos_lat * math.sin(lonr)
    Unit_vector[2] = math.sin(latr)

##############################################################################

# These are all of the major temperature conversions.
def c_to_k(cel):
    return (cel + 273.15)

def k_to_c(kel):
    return (kel - 273.15)

def c_to_f(cel):
    return (cel*1.8+32.0)

def f_to_c(far):
    return ((far-32.0)*5.0/9.0)

def k_to_f(kel):
    return (c_to_f(k_to_c(kel)))

def f_to_k(far):
    return (c_to_k(f_to_c(far)))

###########################################################################

# mmdd_to_jd converts a regular date (YYYYMMDD) to a Julian date (YYYYJdd)
def mmdd_to_jd(date):

    t = (int(date[:4]), int(date[4:6]), int(date[6:8]), 0, 0, 0, 0, 0, -1)
    utime = time.mktime(t)
    t = time.localtime(utime)
    newdate = '%d%03d' % (t[0], t[7])
    return newdate

##########################################################################
    
# jd_to_mmdd converts a Julian date (YYYYJdd) to a regular date (YYYYMMDD)
#  Using the time module gives errors when this is attempted, so it must
#  be forced.
def jd_to_mmdd(jdate):
    
    year = int(jdate[:4])
    jd = int(jdate[4:])
    if ((year%4) == 0):
        if (jd <= 31):
            month = 1
            day = jd
        elif (jd <= 60):
            month = 2
            day = jd - 31
        elif (jd <= 91):
            month = 3
            day = jd - 60
        elif (jd <= 121):
            month = 4
            day = jd - 91
        elif (jd <= 152):
            month = 5
            day = jd - 121
        elif (jd <= 182):
            month = 6
            day = jd - 152
        elif (jd <= 213):
            month = 7
            day = jd - 182
        elif (jd <= 244):
            month = 8
            day = jd - 213
        elif (jd <= 274):
            month = 9
            day = jd - 244
        elif (jd <= 305):
            month = 10
            day = jd - 274
        elif (jd <= 335):
            month = 11
            day = jd - 305
        elif (jd <= 366):
            month = 12
            day = jd - 335
        else:
            print "Not a valid Julian day"
            sys.exit(2)

    elif ((year%4) != 0):
        if (jd <= 31):
            month = 1
            day = jd
        elif (jd <= 59):
            month = 2
            day = jd - 31
        elif (jd <= 90):
            month = 3
            day = jd - 59
        elif (jd <= 120):
            month = 4
            day = jd - 90
        elif (jd <= 151):
            month = 5
            day = jd - 120
        elif (jd <= 181):
            month = 6
            day = jd - 151
        elif (jd <= 212):
            month = 7
            day = jd - 181
        elif (jd <= 243):
            month = 8
            day = jd - 212
        elif (jd <= 273):
            month = 9
            day = jd - 243
        elif (jd <= 304):
            month = 10
            day = jd - 273
        elif (jd <= 334):
            month = 11
            day = jd - 304
        elif (jd <= 365):
            month = 12
            day = jd - 334
        else:
            print "Not a valid Julian Day"
            sys.exit(2)

    newdate = "%d%02d%02d" % (year, month, day)
    return newdate

###########################################################################

# calc_heading calculates a heading given a latitude and longitude pair
def calc_heading(lat1, lon1, lat2, lon2):

    # The equation is set up such that western hemisphere longitudes are positive
    if lon1 < 0.0:
        lon1 = lon1 * -1.0
    if lon2 < 0.0:
        lon2 = lon2 * -1.0

    # Convert lats and lons to radians for the calculation
    lat1 = lat1 * (pi/180.0)
    lon1 = lon1 * (pi/180.0)
    lat2 = lat2 * (pi/180.0)
    lon2 = lon2 * (pi/180.0)

    # Do the calculation.  This comes from http://williams.best.vwh.net/avform.htm and
    # http://mathforum.org/library/drmath/view/55417.html
    num = math.atan2(math.sin(lon1-lon2)*math.cos(lat2), \
                     math.cos(lat1)*math.sin(lat2)-math.sin(lat1)*math.cos(lat2)*math.cos(lon1-lon2))
    tc1 = num % (2*math.pi)

    return (tc1/(math.pi/180.0))

###############################################################################

# dir2num takes a direction (N, NNE, etc.) and returns a compass heading (0-360)
def dir2num(d):

    dir = string.upper(d)
    if dir == 'N':
        return 0.0
    elif dir == 'NNE':
        return 22.5
    elif dir == 'NE':
        return 45.0
    elif dir == 'ENE':
        return 67.5
    elif dir == 'E':
        return 90.0
    elif dir == 'ESE':
        return 112.5
    elif dir == 'SE':
        return 135
    elif dir == 'SSE':
        return 157.5
    elif dir == 'S':
        return 180.0
    elif dir == 'SSW':
        return 202.5
    elif dir == 'SW':
        return 225.0
    elif dir == 'WSW':
        return 247.5
    elif dir == 'W':
        return 270.0
    elif dir == 'WNW':
        return 292.5
    elif dir == 'NW':
        return 315.0
    elif dir == 'NNW':
        return 337.5
    else:
        return -999.9
