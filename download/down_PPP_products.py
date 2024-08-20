#! python3
"""
@Description 	:   PPP-RTK processing/plotting library for POS/TRO/IONO
@Author      	:   Lewen Zhao
@Version     	:   1.0
@Contact     	:   lwzhao@nuist.edu.cn
@Contact     	:   navsense_support@163.com
@OnlineService :	http://1.13.180.60:8800/login2
@Time        	:   2024/03/29 12:13:45
"""
import sys
from datetime import datetime
from dateutil.relativedelta import relativedelta
from down_tools import *
from down_eph_clk import wget_eph_clk
from down_rinex_DSB import wget_dsb
from down_rinex_nav_3x import  wget_rinex3
from down_rinex_nav_4x import  wget_rinex4
        
def down_PPP_data(tstart: datetime, tlen: int, prd_name:str, dir_dst: str):
    tmn = tstart
    if not os.path.exists(dir_dst):
        os.makedirs(dir_dst)
    for day in range(tlen):
        wget_eph_clk(tmn, dir_dst, prd_name)
        wget_rinex3(tmn, dir_dst)
        wget_rinex4(tmn, dir_dst)
        wget_dsb(tmn, dir_dst)
        tmn = tmn + relativedelta(days=1)
        
if __name__ == "__main__":
    tstart = datetime(2024, 8, 16, 0, 0, 0)
    tlen = 1
    dir_dst = r'D:\dataset_PPP\test_ppp_rtk_service'
    down_PPP_data(tstart, tlen,"WUM", dir_dst)
