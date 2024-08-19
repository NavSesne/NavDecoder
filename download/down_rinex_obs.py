"""
@Description 	:   PPP-RTK processing/plotting library for POS/TRO/IONO
@Author      	:   Lewen Zhao
@Version     	:   1.0
@Contact     	:   lwzhao@nuist.edu.cn
@Contact     	:   navsense_support@163.com
@OnlineService :	http://1.13.180.60:8800/login2
@Time        	:   2024/03/29 12:14:20
"""

import datetime
import os
from dateutil.relativedelta import relativedelta
from down_tools import *
from concurrent.futures import ThreadPoolExecutor
import timeit

def getLongRinexName(station):
    mapping_file=r'D:\lewen\source_code\PPPRTK_src\bin\RINEX_9CHAR'
    matches = ""
    with open(mapping_file, 'r') as file:
        for line in file:
            if line[:4] == station:
                matches=line.strip()
            
    return matches

def getRemoteURL(tmn,ftp_name):
    yyyy=tmn.year
    doy  = ymd2doy(tmn.year, tmn.month, tmn.day)
    yy = tmn.year-2000 if tmn.year>=2000 else tmn.year-1900 


    URL_CDDIS = f"ftps://gdc.cddis.eosdis.nasa.gov/pub/gnss/data/daily/{yyyy}/{doy:03d}/{yy}d/"
    URL_IGN = f"ftp://igs.ign.fr/pub/igs/data/{yyyy}/{doy:03d}/"
    URL_WHU = f"ftp://igs.gnsswhu.cn/pub/gps/data/daily/{yyyy}/{doy:03d}/{yy}d/"

    rpath=""
    if ftp_name=="WHU":
        rpath=URL_WHU
    if ftp_name=="CDDIS":
        rpath=URL_CDDIS
    if ftp_name=="IGN":
        rpath=URL_IGN
    
    return rpath
    
def down_rinex_obs(tmn,local_dir,station,ftp_name):

    yyyy=tmn.year
    doy  = ymd2doy(tmn.year, tmn.month, tmn.day)
    rpath=getRemoteURL(tmn,ftp_name)

    name9=getLongRinexName(station)
    name = f'{name9}_R_{yyyy}{doy:03d}0000_01D_30S_MO'
    local_rnx=os.path.join(local_dir,name+".rnx")
    if not os.path.exists(local_rnx):
        local_file=download_compress(rpath,local_dir,name+'.crx','.gz')
        cmd_command = f'{bin_crx}  {local_file}'
        os.system(cmd_command)
        if os.path.exists(local_file):
            os.remove(local_file)

def down_rinex_obs_highrate(tmn,local_dir,station,ftp_name):
    yyyy=tmn.year
    doy  = ymd2doy(tmn.year, tmn.month, tmn.day)
    rpath=getRemoteURL(tmn,ftp_name)
    yy = tmn.year-2000 if tmn.year>=2000 else tmn.year-1900 

    URL_CDDIS = f"ftps://gdc.cddis.eosdis.nasa.gov/pub/gnss/data/highrate/{yyyy}/{doy:03}/{yy}d/"
    URL_IGN = f"ftp://igs.ign.fr/pub/igs/data/highrate/{yyyy}/{doy:03}/"
    URL_WHU = f"ftp://igs.gnsswhu.cn/pub/highrate/{yyyy}/{doy:03}/{yy}d/"

    rpath=""
    if ftp_name=="WHU":
        rpath=URL_WHU
    if ftp_name=="CDDIS":
        rpath=URL_CDDIS
    if ftp_name=="IGN":
        rpath=URL_IGN

    name9=getLongRinexName(station)
    # generate the hourly name
    hourly_files=[]
    for hour in range(24):
        rpath_hour=f"{rpath}{hour:02}/"
        for min in [0,15,30,45]:
            name = f'{name9}_R_{yyyy}{doy:03d}{hour:02}{min:02}_15M_01S_MO'
            local_rnx=os.path.join(local_dir,name+".rnx")
            hourly_files.append(local_rnx)
            if not os.path.exists(local_rnx):
                local_file=download_compress(rpath_hour,local_dir,name+'.crx','.gz')
                cmd_command = f'{bin_crx}  {local_file}'
                os.system(cmd_command)
                if os.path.exists(local_file):
                    os.remove(local_file)
    merged_files=f'{name9}_R_{yyyy}{doy:03d}????_15M_01S_MO.rnx'
    name_1s=f'{name9}_R_{yyyy}{doy:03d}0000_01D_01S_MO.rnx'
    name_5s=f'{name9}_R_{yyyy}{doy:03d}0000_01D_05S_MO.rnx'
    cmd_merge="%s -finp %s -fout %s -kv -splice_memsave -satsys GEC" %(bin_gfzrnx,os.path.join(local_dir,merged_files),os.path.join(local_dir,name_1s))
    cmd_sample="%s -finp %s -fout %s -smp 5 -satsys GEC" %(bin_gfzrnx,os.path.join(local_dir,name_1s),os.path.join(local_dir,name_5s))

    os.system(cmd_merge)
    os.system(cmd_sample)

    # delete the houly files
    for file in hourly_files:
        if os.path.exists(file):
            os.remove(file)


if __name__ == "__main__":
    tstart = datetime.datetime(2024, 6,  6, 0, 0, 0)
    tlen=30
    # stations=['HKSL','HKWS','NCKU','KMNM','TNML','TWTF','LHAZ','JFNG','WUHN','SHAO','BJFS','BJNM','CHAN','URUM']
    stations=["WUH2","GAMG","ULAB","URUM","POL2","BIK0","LCK4","JFNG","WUH2"] #JFNG
    data_center="WHU"  #WHU  IGN CDDIS
    local_dir = r'E:\GNSS_Data\MGEX_OBS'
    # local_dir = r'E:\GNSS_Data\MGEX_OBS\highrate'
    dates=[]
    tmn=tstart
    with ThreadPoolExecutor(max_workers=5) as executor:
        start_time = timeit.default_timer()
        for day in range(tlen):
            futures = []
            for station in stations:
                future = executor.submit(down_rinex_obs_highrate, tmn, local_dir, station, data_center)
                # futures.append(future)
                # data center: WHU CDDIS
                # down_rinex_obs(tmn,local_dir,station,data_center)
                # down_rinex_obs_highrate(tmn,local_dir,station,data_center)
            tmn = tmn + relativedelta(days=1)
        end_time = timeit.default_timer() - start_time
        print(f"============下载第{tmn.strftime('%Y-%m-%d %H:%M:%S')}天数据完成! 总共耗时: {end_time:.2f} 秒")
