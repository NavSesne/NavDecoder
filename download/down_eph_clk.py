#! python3
#coding:utf-8
"""
@Description 	:   PPP-RTK processing/plotting library for POS/TRO/IONO
@Author      	:   Lewen Zhao
@Version     	:   1.0
@Contact     	:   lwzhao@nuist.edu.cn
@Contact     	:   navsense_support@163.com
@OnlineService :	http://1.13.180.60:8800/login2
@Time        	:   2024/03/29 12:13:18
"""

import datetime,os
from dateutil.relativedelta import relativedelta
from down_tools import *
from download.cmn_tools import *
from concurrent.futures import ThreadPoolExecutor, as_completed
URL_IGN = "ftp://igs.ign.fr/pub/igs/products/mgex/"
URL_WHU = "ftp://igs.gnsswhu.cn/pub/gnss/products/mgex/"
URL_Default = URL_WHU
URL_CNT = "http://www.ppp-wizard.net/products/REAL_TIME/"
prd_enum_set = {"IGS", "GRM", "COM", "WUM", "GFR", "WHR", "CNT"}
def generate_filenames_and_url(center, tmn):
    # Check if the center is valid
    if center not in prd_enum_set:
        raise ValueError(f"Center '{center}' is not recognized. Please use one of the following: {', '.join(prd_enum_set)}")
    # 根据center和时间t生成文件名和下载URL
    # 此处需要根据实际情况填充完整的逻辑
    year= tmn.year  # datetime 对象的年份属性
    time0=epoch2time([tmn.year, tmn.month, tmn.day, tmn.hour, tmn.minute, tmn.second])
    week,dow=time2gpst(time0)
    doy=time2doy(time0)
    ephf_name = clkf_name = biaf_name = suffix=""
    download_url = ""
    current_date = datetime.now()
    ddays = (current_date - tmn).days
    dataNRT,dataRAP,dataFIN=False,False,False
    default=None
    if ddays<1:
        dataNRT=True
    elif 1<ddays and ddays<7:  #Rapied products
        dataNRT=True
        dataRAP=True
    elif ddays>=7: #Final products
        dataNRT=True
        dataRAP=True
        dataFIN=True
    if dataNRT or dataRAP or dataFIN: #Real-time prodcuts from B2b and HAS
        if center == "B2B" or center == "HAS":
            print("Download real-time archive of PPP-B2b and Galileo HAS products")
        elif center=="WUM_NRT" or (dataNRT==True and dataRAP==False and dataFIN==False):
            t_NRT=tmn+relativedelta(days=-1) # the orbit and clock in the file starts one-day ago
            time0=epoch2time([t_NRT.year, t_NRT.month, t_NRT.day, t_NRT.hour, t_NRT.minute, t_NRT.second])
            week0,dow0=time2gpst(time0)
            doy0=time2doy(time0)
            ephf_name = f"WUM0MGXNRT_{year}{doy0:03d}{t_NRT.hour:02d}00_02D_05M_ORB.SP3"
            clkf_name = f"WUM0MGXNRT_{year}{doy0:03d}{t_NRT.hour:02d}00_02D_30S_CLK.CLK"
            biaf_name = f"WUM0MGXNRT_{year}{doy0:03d}0000_01D_01D_OSB.BIA"
            download_url = f"{URL_Default}{week0}/"
            suffix='.gz'
            if (dataNRT==True and dataRAP==False and dataFIN==False):
                print(f"Specificted products not found, {center}==>WUM_NRT")

    if dataRAP or dataFIN:  #Rapied products
        if center == "GFR":
            ephf_name = f"GFZ0MGXRAP_{year}{doy:03d}0000_01D_05M_ORB.SP3"
            clkf_name = f"GFZ0MGXRAP_{year}{doy:03d}0000_01D_30S_CLK.CLK"
            download_url = f"{URL_Default}{week}/"
            suffix='.gz'
        elif center == "CNT":
            ephf_name = f"cnt{week}{dow}.sp3"
            clkf_name = f"cnt{week}{dow}.clk"
            biaf_name = f"cnt{week}{dow}.bia"
            URL_CNT ="http://www.ppp-wizard.net/products/REAL_TIME/"
            download_url = URL_CNT
            suffix='.gz'        
        elif center == "WHR" or (dataRAP==True and dataFIN==False):
            ephf_name = f"WUM0MGXRAP_{year}{doy:03d}0000_01D_05M_ORB.SP3"
            clkf_name = f"WUM0MGXRAP_{year}{doy:03d}0000_01D_30S_CLK.CLK"
            biaf_name = f"WUM0MGXRAP_{year}{doy:03d}0000_01D_01D_OSB.BIA"
            URL_WHU="ftp://igs.gnsswhu.cn/pub/whu/phasebias/"
            download_url = f"{URL_WHU}{year}/"
            suffix='.gz'
            if (dataRAP==True and dataFIN==False):
                print(f"Specificted products not found, {center}==>WHR")

    if dataFIN: #Final products
        if center == "GRM":
            ephf_name = f"GRG0MGXFIN_{year}{doy:03d}0000_01D_05M_ORB.SP3"
            clkf_name = f"GRG0MGXFIN_{year}{doy:03d}0000_01D_30S_CLK.CLK"
            download_url = f"{URL_Default}{week}/"
            suffix='.gz'
        elif center == "IGS":
            ephf_name = f"igs{year}{doy:03d}.sp3"
            clkf_name = f"igs{year}{doy:03d}.clk_30s"
            URL_IGS ="ftp://igs.gnsswhu.cn/pub/gps/products/"
            download_url = f"{URL_IGS}{week}/"
            suffix='.Z'
        elif center == "WUM" :
            ephf_name = f"WUM0MGXFIN_{year}{doy:03d}0000_01D_05M_ORB.SP3"
            clkf_name = f"WUM0MGXFIN_{year}{doy:03d}0000_01D_30S_CLK.CLK"
            download_url = f"{URL_Default}{week}/"
            suffix='.gz'

    return ephf_name, clkf_name, biaf_name, download_url,suffix

def wget_eph_clk(tmn,local_dir,center):
    doy  = ymd2doy(tmn.year, tmn.month, tmn.day)
    yy = tmn.year-2000 if tmn.year>=2000 else tmn.year-1900 
    name = 'BRD400DLR_S_%4d%03d0000_01D_MN.rnx' % (tmn.year, doy)
    compressed_name = name + '.gz'
    rpath = 'ftp://igs.gnsswhu.cn/pub/gps/data/daily/%d/%03d/%2dp/' % (tmn.year, doy, tmn.year % 100)
    ephf_name, clkf_name, biaf_name, remote_path,suffix = generate_filenames_and_url(center, tmn)

    if not os.path.exists(local_dir):
        os.makedirs(local_dir)

    if ephf_name:
        download_url=remote_path
        if "whu/phasebias"in download_url:
            download_url = f"{download_url}orbit/"
        download_compress(download_url,local_dir, ephf_name,suffix)
    if clkf_name:
        download_url=remote_path
        if "whu/phasebias"in download_url:
            download_url = f"{download_url}clock/"
        download_compress(download_url,local_dir, clkf_name,suffix)
    if biaf_name:
        download_url=remote_path
        if "whu/phasebias"in download_url:
            download_url = f"{download_url}bias/"
        download_compress(download_url,local_dir, biaf_name,suffix)

if __name__ == "__main__":
    tstart = datetime(2024, 5, 5, 0, 0, 0)
    tlen=60
    dir_dst = r'E:\GNSS_Data\products\eph'
    tmn=tstart

    max_workers = 5  # 根据需要调整最大工作线程数
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for day in range(tlen):
            future = executor.submit(wget_eph_clk, tmn, dir_dst, "WUM")
            futures.append(future)
            tmn = tmn + relativedelta(days=1)

        # 等待所有任务完成
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"下载任务发生错误: {e}")

    print("所有下载任务完成!")