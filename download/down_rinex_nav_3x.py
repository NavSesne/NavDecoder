#! python3
"""
@Description 	:   PPP-RTK processing/plotting library for POS/TRO/IONO
@Author      	:   Lewen Zhao
@Version     	:   1.0
@Contact     	:   navsense_support@163.com
@OnlineService :	http://1.13.180.60:8800/login2
@Time        	:   2024/03/20 23:47:21
"""

import datetime
from dateutil.relativedelta import relativedelta
from down_tools import *

def wget_rinex3(tmn,local_dir):
    doy  = ymd2doy(tmn.year, tmn.month, tmn.day)
    yy = tmn.year-2000 if tmn.year>=2000 else tmn.year-1900 
    name = 'BRDC00GOP_R_%4d%03d0000_01D_MN.rnx' % (tmn.year, doy)
    rpath = 'ftp://ftp.pecny.cz/LDC/orbits_brd/gop3/%d/' % (tmn.year)
    download_compress(rpath,local_dir,name,'.gz')

        
if __name__ == "__main__":
    tstart = datetime.datetime(2024,  4,  8, 0, 0, 0)
    tlen=2
    dir_dst = r'D:\2024-04-07-debug'
    tmn=tstart
    for day in range(tlen):
        wget_rinex3(tmn,dir_dst)
        tmn = tmn + relativedelta(days=1)
