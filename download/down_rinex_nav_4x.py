#! python3
"""
@Description 	:   PPP-RTK processing/plotting library for POS/TRO/IONO
@Author      	:   Lewen Zhao
@Version     	:   1.0
@Contact     	:   lwzhao@nuist.edu.cn
@Contact     	:   navsense_support@163.com
@OnlineService :	http://1.13.180.60:8800/login2
@Time        	:   2024/03/29 12:13:06
"""

import datetime
from dateutil.relativedelta import relativedelta
from down_tools import *

def wget_rinex4(tmn,local_dir):
    doy  = ymd2doy(tmn.year, tmn.month, tmn.day)
    yy = tmn.year-2000 if tmn.year>=2000 else tmn.year-1900 
    name = 'BRD400DLR_S_%4d%03d0000_01D_MN.rnx' % (tmn.year, doy)
    rpath = 'ftp://igs.gnsswhu.cn/pub/gps/data/daily/%d/%03d/%2dp/' % (tmn.year, doy, tmn.year % 100)
    download_compress(rpath,local_dir,name,'.gz')

        
if __name__ == "__main__":
    tstart = datetime.datetime(2024,  3,  19, 0, 0, 0)
    tlen=1
    dir_dst = r'E:\GNSS_Data\B2B_CORR\product'
    tmn=tstart
    for day in range(tlen):
        wget_rinex4(tmn,dir_dst)
        tmn = tmn + relativedelta(days=1)
