
"""
@Description 	:   PPP-RTK processing/plotting library for POS/TRO/IONO
@Author      	:   Lewen Zhao
@Version     	:   1.0
@Contact     	:   lwzhao@nuist.edu.cn
@Contact     	:   navsense_support@163.com
@OnlineService :	http://1.13.180.60:8800/login2
@Time        	:   2024/03/29 12:14:34
"""
#coding=utf-8
import os
import sys
from rtkcmn import *
import datetime
from dateutil.relativedelta import relativedelta
from down_tools import *

def get_refCoor(site,rcvRef,xyzRef):
    ref=[]
    if len(site)>4:
        site=site[0:4]
    for ii in range(len(rcvRef)):
        sites=rcvRef[ii]
        if(sites.lower()==site.lower()):
            ref=[float(xyzRef[ii][0]),float(xyzRef[ii][1]),float(xyzRef[ii][2])]

    return ref

#if __name__ == "__main__":
def get_xyz_Sinex(year,month,day,acs,sites,dir_dest):
    xyzRef=[]
    rcvRef=[]
    rpath,fname,suffix="","",""
    try: 
        #year,month,day,acs,sites=sys.argv[1:6]
        time0=epoch2time([year,month,day,0,0,0])
        week,dow=time2gpst(time0)
        doy=time2doy(time0)
        #URL_CDDIS="https://cddis.nasa.gov/archive/gnss/products/2290/"
        URL_IGN="ftp://igs.ensg.ign.fr/pub/igs/products"
        URL_WHU="ftp://igs.gnsswhu.cn/pub/gps/products/"
        if acs.upper()=="IGS":
            #fname='igs%02dP%04d%01d.snx' %(int(year)-2000,week,dow)
            fname='igs%02dP%04d.snx' %(int(year)-2000,week)
            suffix='.Z'
            rpath='%s/%04d/' %(URL_IGN,week)
        elif  acs.upper()=="COD"  :
            #fname='cod%04d%01d.snx' %(week,dow)
            fname='igs%02dP%04d.snx' %(int(year)-2000,week)
            suffix='.Z'
            rpath='%s/%04d/' %(URL_IGN,week)
        elif  acs.upper()=="WHU"  :
            fname='IGS0OPSSNX_%4d%03d0000_01D_01D_SOL.SNX' %(year,doy)
            suffix='.gz'
            rpath='%s/%04d/' %(URL_WHU,week)
        else:
            return 

        if len(dir_dest)==0:
            dir_dest=sys.path[0]

        snx_file=download_compress(rpath,dir_dest,fname,suffix)

        with open(snx_file, 'r', encoding='utf-8') as fo:
            lines = fo.readlines()
        
        site_info = {site.upper(): [] for site in sites}  # 使用字典存储每个站点的信息

        for line in lines:
            site_code=None
            parts = line.strip().split()
            if not parts:
                continue  # 跳过空行
            if parts[0] in site_info:
                site_code = parts[0]
            elif len(parts) > 2 and parts[2] in site_info:
                site_code = parts[2]

            # 如果找到站点代码，则处理该行
            if site_code:
                site_info[site_code].append(line)
        for ss, staInfo in site_info.items():
            if not staInfo:
                print(f"Stations not in the sinex {ss}")
                continue
            if len(staInfo) != 11:
                print(f"Format not standard for station {ss}")
                print("staInfo ", staInfo)
                continue
            
            rcvt = staInfo[1][42:62].strip()
            antt = staInfo[2][42:62].strip()
            strUNE = staInfo[3].split()
            dneu = [float(strUNE[-2]), float(strUNE[-1]), float(strUNE[-3])]
            dxyz = [float(staInfo[-3].split()[8]), float(staInfo[-2].split()[8]), float(staInfo[-1].split()[8])]

            xyzRef.append(dxyz)
            rcvRef.append(ss)
            blh = xyz2blh(*dxyz)  # 假设xyz2blh已经定义
            # print(f"{ss}  {blh[0]} {blh[1]} {blh[2]}")
    except Exception as e:
        print(sys.argv)
        print(e)
    return rcvRef,xyzRef

if __name__ == "__main__":
    tstart = datetime.datetime(2023,  11,  20, 0, 0, 0)
    tlen=11
    dir_dst = r'D:\lewen\Work_NUIST\RTPPP_Eval\product'
    tmn=tstart
    sites=["JFNG"]
    for day in range(tlen):
        rcvRef, xyzRef = get_xyz_Sinex(tstart.year, tstart.month, tstart.day, "WHU", sites, dir_dst)
        tmn = tmn + relativedelta(days=1)
