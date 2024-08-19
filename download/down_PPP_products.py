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

from datetime import timedelta
from dateutil.relativedelta import relativedelta
from down_tools import *
from down_eph_clk import wget_eph_clk
from down_rinex_DSB import wget_dsb
from down_rinex_nav_3x import  wget_rinex3
from down_rinex_nav_4x import  wget_rinex4
        
def year_doy_to_date(year, doy):
    # Ensure year is an integer
    year = int(year)
    # Ensure doy (day of year) is an integer
    doy = int(doy)
    # Create a datetime object for the first day of the year
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    # Return the date in 'YYYY-MM-DD' format
    return date


def down_PPP_data_doy(year,doy,tlen,dir_dst):
    tmn = year_doy_to_date(year,doy)
    if not os.path.exists(dir_dst):
        os.makedirs(dir_dst)
    for day in range(tlen):
        wget_eph_clk(tmn, dir_dst, "WHR")
        wget_rinex3(tmn, dir_dst)
        wget_rinex4(tmn, dir_dst)
        wget_dsb(tmn, dir_dst)
        tmn = tmn + relativedelta(days=1)

def down_PPP_data_ymd(tstart: datetime, tlen: int, dir_dst: str):
    tmn = tstart
    if not os.path.exists(dir_dst):
        os.makedirs(dir_dst)
    for day in range(tlen):
        wget_eph_clk(tmn, dir_dst, "WHR")
        wget_rinex3(tmn, dir_dst)
        wget_rinex4(tmn, dir_dst)
        wget_dsb(tmn, dir_dst)
        tmn = tmn + relativedelta(days=1)
        
if __name__ == "__main__":
    # down_PPP_data_ymd(datetime(2024, 5, 1, 0, 0, 0)  , 60, r'E:\GNSS_Data\products\eph')
    down_PPP_data_doy(2024,134,2,r"D:\work_lewen\source_code\git_lewen\NavDecoder\test_data")