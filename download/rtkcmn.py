"""
@Description   : PPP-RTK processing/plotting library for POS/TRO/IONO
@Author        : Lewen Zhao
@Version       : 1.0
@Contact       : lwzhao@nuist.edu.cn
@Contact       : navsense_support@163.com
@OnlineService : http://1.13.180.60:8800/login2
@Time :   2024/07/17 13:02:20
"""
#!/usr/bin/env python3
#coding=utf-8
import math
from copy import copy, deepcopy
from math import floor, sin, cos, sqrt, asin, atan2, fabs
from datetime import datetime
A = 6378137.0
B = 6356752.3142
gpst0 = [1980, 1, 6, 0, 0, 0]

def leaps(tgps):
    """ return leap seconds (TBD) """
    return -18.0

class gtime_t():
    """ class to define the time """

    def __init__(self, time=0, sec=0.0):
        self.time = time
        self.sec = sec

def xyz2blh(x, y, z):
    """Convert XYZ coordinates to BLH,
    return tuple(latitude, longitude, height).
    """
    e = math.sqrt(1 - (B**2)/(A**2))
    # calculate longitude, in radians
    if x >= 0:
        longitude = math.atan(y/x)
    elif y >= 0:
        longitude = math.pi + math.atan(y/x)
    else:
        longitude = -math.pi + math.atan(y/x)
    # calculate latitude, in radians
    xy_hypot = math.hypot(x, y)

    lat0 = 0
    latitude = math.atan(z / xy_hypot)

    while abs(latitude - lat0) > 1E-5:
        lat0 = latitude
        N = A / math.sqrt(1 - e**2 * math.sin(lat0)**2)
        latitude = math.atan((z + e**2 * N * math.sin(lat0)) / xy_hypot)
    # calculate height, in meters
    height = z / math.sin(latitude) - N * (1 - e**2)
    # convert angle unit to degrees
    longitude = math.degrees(longitude)
    latitude = math.degrees(latitude)

    return latitude, longitude, height


def xyz2neu(x0, y0, z0, x, y, z):
    """
    Convert cartesian coordinate system to site-center system.
    Input paraments:
    - x0, y0, z0: coordinate of centra site,
    - x, y, z: coordinate to be converted.
    """
    # calculate the lat, lon and height of center site
    lat, lon, _ = xyz2blh(x0, y0, z0)
    # convert angle unit to radians
    lat, lon = math.radians(lat), math.radians(lon)
    # calculate NEU
    north = (-math.sin(lat) * math.cos(lon) * (x - x0) -
             math.sin(lat) * math.sin(lon) * (y - y0) +
             math.cos(lat) * (z - z0))
    east = -math.sin(lon) * (x - x0) + math.cos(lon) * (y - y0)
    up = (math.cos(lat) * math.cos(lon) * (x- x0) +
          math.cos(lat) * math.sin(lon) * (y - y0) +
          math.sin(lat) * (z - z0))

    return north, east, up

def blh2xyz(latitude, longitude, height):
    """Convert BLH coordinates to XYZ.
    return tuple(X, Y, Z).
    """
    # convert angle unit to radians
    latitude = math.radians(latitude)
    longitude = math.radians(longitude)

    e = math.sqrt(1 - (B**2)/(A**2))
    N = A / math.sqrt(1 - e**2 * math.sin(latitude)**2)
    # calculate X, Y, Z
    X = (N + height) * math.cos(latitude) * math.cos(longitude)
    Y = (N + height) * math.cos(latitude) * math.sin(longitude)
    Z = (N * (1 - e**2) + height) * math.sin(latitude)

    return X, Y, Z

def epoch2time(ep):
    """ calculate time from epoch """
    doy = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
    time = gtime_t()
    year = int(ep[0])
    mon = int(ep[1])
    day = int(ep[2])

    if year < 1970 or year > 2099 or mon < 1 or mon > 12:
        return time
    days = (year-1970)*365+(year-1969)//4+doy[mon-1]+day-2
    if year % 4 == 0 and mon >= 3:
        days += 1
    sec = int(ep[5])
    time.time = days*86400+int(ep[3])*3600+int(ep[4])*60+sec
    time.sec = ep[5]-sec
    return time

def gpst2utc(tgps, leaps_=-18):
    """ calculate UTC-time from gps-time """
    tutc = timeadd(tgps, leaps_)
    return tutc

def utc2gpst(tutc, leaps_=-18):
    """ calculate UTC-time from gps-time """
    tgps = timeadd(tutc, -leaps_)
    return tgps

def timeadd(t: gtime_t, sec: float):
    """ return time added with sec """
    tr = copy(t)
    tr.sec += sec
    tt = floor(tr.sec)
    tr.time += int(tt)
    tr.sec -= tt
    return tr

def timediff(t1: gtime_t, t2: gtime_t):
    """ return time difference """
    dt = t1.time - t2.time
    dt += (t1.sec - t2.sec)
    return dt

def gpst2time(week, tow):
    """ convert to time from gps-time """
    t = epoch2time(gpst0)
    if tow < -1e9 or tow > 1e9:
        tow = 0.0
    t.time += 86400*7*week+int(tow)
    t.sec = tow-int(tow)
    return t

def gpst2datetime(week, tow):
    """ convert to time from python date-time """
    t = epoch2time(gpst0)
    if tow < -1e9 or tow > 1e9:
        tow = 0.0
    t.time += 86400*7*week+int(tow)
    t.sec = tow-int(tow)
    ep=time2epoch(t)
    ep[5]=int(ep[5])
    return datetime(*ep)

def datetime2time(tmn):
    time0=epoch2time([tmn.year, tmn.month, tmn.day, tmn.hour, tmn.minute, tmn.second])
    return time0

def time2datetime(inp_time):
    ep=time2epoch(inp_time)

    return datetime(*ep)

def time2gpst(t: gtime_t):
    """ convert to gps-time from time """
    t0 = epoch2time(gpst0)
    sec = t.time-t0.time
    week = int(sec/(86400*7))
    tow = sec-week*86400*7+t.sec
    return week, tow

def time2epoch(t):
    """ convert time to epoch """
    mday = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31,
            30, 31, 31, 30, 31, 30, 31, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31,
            30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    days = int(t.time/86400)
    sec = int(t.time-days*86400)
    day = days % 1461
    for mon in range(48):
        if day >= mday[mon]:
            day -= mday[mon]
        else:
            break
    ep = [0, 0, 0, 0, 0, 0]
    ep[0] = 1970+days//1461*4+mon//12
    ep[1] = mon % 12+1
    ep[2] = day+1
    ep[3] = sec//3600
    ep[4] = sec % 3600//60
    ep[5] = sec % 60+t.sec
    return ep

def time2doy(t):
    """ convert time to epoch """
    ep = time2epoch(t)
    ep[1] = ep[2] = 1.0
    ep[3] = ep[4] = ep[5] = 0.0
    return int(timediff(t, epoch2time(ep))/86400+1)
