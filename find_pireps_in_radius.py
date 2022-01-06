#!/usr/local/python/bin/python

# This program takes in a lat/lon (in decimal degrees) and finds
# all of the PIREPs that are within a given radius of the point
# for a certain start date and number of dates and outputs them
# to the output file 

import os, string, sys, getopt, math, time
from Numeric import *
from conv import *

def usage():
    print "usage: %s lat lon radius start_date num_days outfile" % sys.argv[0]
    print "lat and lon are in decimal degrees and represent the center point"
    print "radius is the max distance (in km) that the PIREPs can be from the point"
    print "start_date is the date to start looking for PIREPs (YYYYMMDD)"
    print "num_days is the number of days to look over"
    print "outfile is where the matching PIREPs should go"
    sys.exit(2)

arg_len = len(sys.argv)
if arg_len != 7:
   usage()


clat = float(sys.argv[1])
clon = float(sys.argv[2])
radius = int(sys.argv[3])
start = sys.argv[4]
num_days = int(sys.argv[5])
outfile = sys.argv[6]

if clon < 0.0:
    clon = clon * -1.0
    
# Open the outfile
fout = open(outfile, "w")

t = (int(start[:4]), int(start[4:6]), int(start[6:]), 12, 0, 0, 0, 0, 0)
itime = time.mktime(t)

n = 0
for i in xrange(0,num_days):
   # Convert the new Unix time into a local time to get the next day
   utime = itime + (i * 86400)  # 86400 seconds per day
   t = time.localtime(utime)
   year = t[0]
   ymd = "%d%02d%02d" % (year, t[1], t[2])
   print "Date being processed: %s\n" % ymd

   # Open each PIREP file and go through it
   for h in xrange(0,23):
       pfile = "/var/autofs/mnt/khufu_d1/cwolff/decoded_pireps/%s%02d_pireps" % (ymd, h)
       try:
           pin = open(pfile, "r")
       except:
           try:
               pfile = "/var/autofs/mnt/nikara_d3/cwolff/pireps/decoded_pireps/%s%02d_pireps" % (ymd, h)
               pin = open(pfile, "r")
           except:
               continue

       plines = pin.readlines()

       U = []
       Lat = []
       Lon = []
       T = []
       for line in plines:
           p = string.split(line)
           u = int(p[0])
           plat = float(p[1])
           plon = float(p[2])
           type = p[33]
           ice = int(p[17])

           if ice < -1:
               continue
           
           # Do some error checking
           if plat < 15.0:
               continue
           if plon < 0.0:
               plon = plon * -1.0
           else:
               continue

           # Check to see if we've seen this PIREP before
           U.append(u)
           Lat.append(plat)
           Lon.append(plon)
           T.append(type)
           found = 0
           l = len(U) - 2
           while l >= 0:
               if ((abs(U[l]-u) < 180) and (Lat[l] == plat) and \
                   (Lon[l] == plon) and (T[l] == type)):
                   found = 1
                   break
               l = l - 1
           if found == 1:
               continue

           # Calculate the distance
           Uv1 = [0.0,0.0,0.0]
           Uv2 = [0.0,0.0,0.0]
           defuvp(clat,clon,Uv1)
           defuvp(plat,plon,Uv2)
           dist = abetuv(Uv1,Uv2)/1000.0 # Convert to km

           # Figure out if it's within the radius
           if dist <= radius:
               fout.write(line)
               n = n + 1

       pin.close()
       
fout.close()
print "\nThere are %d icing PIREPs within %d km of %.2f N, %.2f W between %s and %s \n" % \
      (n, radius, clat, clon, start, ymd)
