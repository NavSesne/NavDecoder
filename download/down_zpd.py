#! python3
"""
@Description 	:   PPP-RTK processing/plotting library for POS/TRO/IONO
@Author      	:   Lewen Zhao
@Version     	:   1.0
@Contact     	:   lwzhao@nuist.edu.cn
@Contact     	:   navsense_support@163.com
@OnlineService :	http://1.13.180.60:8800/login2
@Time        	:   2024/03/29 12:14:48
"""
import datetime
import os, sys, glob
from dateutil.relativedelta import relativedelta

def ymd2doy(year, mon, day):
    dn = datetime.datetime(year, mon, day, 0, 0, 0)
    return int(dn.strftime("%j"))


def call_wget_(dir_dst, url):
    if not os.path.exists(dir_dst):
        os.makedirs(dir_dst)
    cmd = 'wget  -P %s %s' % (dir_dst, url)
    print (cmd)
    os.system(cmd)

def call_curl_(dir_dst, url):
    if not os.path.exists(dir_dst):
        os.makedirs(dir_dst)
    os.chdir(dir_dst)
    cmd = 'curl  -c .cookie -n -L -O %s' % (url)
    print (cmd)
    os.system(cmd)

    
def wget_ZPD(tmn,sites):
    rinexo=''
    for site in sites:
        doy  = ymd2doy(tmn.year, tmn.month, tmn.day)
        yy = tmn.year-2000 if tmn.year>=2000 else tmn.year-1900 
        name  = '%4s%03d0.%02dzpd' % (site.lower(), doy, yy)
        name  = 'IGS0OPSFIN_%4d%03d0000_01D_05M_%s00CHN_TRO.TRO' %(tmn.year,doy,site.upper())
        rpath = 'ftp://igs.gnsswhu.cn/pub/gps/products/troposphere/new/%04d/%03d/' % (tmn.year,doy)
        # rpath = 'https://cddis.nasa.gov/archive/gnss/products/troposphere/zpd/%04d/%03d/' % (tmn.year,doy)
        obsfs = '%s%s.gz' %(rpath,name)
        local='%s%s' %(dir_dst,name)
        locals=glob.glob(local)
        if locals.__len__() > 0:        
            continue

        # call_wget_(dir_dst, obsfs)
        call_curl_(dir_dst, obsfs)        
        # local_gz='%s.gz' %(local)
        # cmd='gzip -d %s' %(local_gz)
        # os.system(cmd)
        
if __name__ == "__main__":
    tstart = datetime.datetime(2023,  1,  1, 0, 0, 0)
    #sites =['ALIC', 'COCO', 'DRAO', 'HOB2', 'KIRU', 'MAJU', 'MAS1', 'MAT1', 'MIZU', 'NICO', 'OUS2', 'PIE1', 'RIO2', 'SUTM', 'THTG', 'ULAB', 'VOIM', 'YAR2', 'YEL2', 'ZIMJ']
    #sites =['AREG','ASCG','BOGT','BRST','BSHM','CAS1','CCJ2','CHPG','CKIS','CORD','CPVG','CUSV','DGAV','DUND','FTNA','GANP','GMSD','GODN','GUAM','HARB','HUEG','JFNG','JOG2','KARR','KERG','KIRI','KOKV','LEIJ','LMMF','MADR','MAO0','MARS','MAYG','METG','MIKL','MIZU','NKLG','','NRMD','OHI2','ONS1','PALV','PERT','PIE1','PNGM','POHN','POLV','PTVL','REUN','RGDG','SEYG','SFER','STFU','STK2','TASH','UNB3','UNSA','WARK']
    #sites =['BJFS','CHAN', 'LHAZ', 'SHAO', 'URUM', 'WUHN']
    sites =['BJFS']
    # sites =["YEL2","REYK","KZN2","ULAB","LHAZ","PNGM","NNOR","KERG","NKLG","CPVG","BIK0","WROC","ALGO","LMMF","POVE","JPLM","FAA1","CHPG","RGDG","BJFS","CHAN","LHAZ","URUM","WUHN","WUSH"]
    tlen=365
    dir_dst = r'D:\Data_lewen\data\ppp_cmonoc\IGS_ZPD_tmp\\'
    tmn=tstart
    for day in range(tlen):
        wget_ZPD(tmn,sites)
        tmn = tmn + relativedelta(days=1)        
