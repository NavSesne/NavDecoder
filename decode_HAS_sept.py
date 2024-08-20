"""
    decoder for Galileo HAS products
"""

import sys, os
from binascii import unhexlify
import bitstruct as bs
import numpy as np
import os
from copy import deepcopy
from tqdm import tqdm
from B2b_HAS_decoder.gnss import *
from B2b_HAS_decoder.peph import peph
from B2b_HAS_decoder.cssr_has_sept import cssr_has
from B2b_HAS_decoder.rinex import rnxdec
from datetime import datetime, timedelta
from B2b_HAS_decoder.cssrlib import sCSSR,sCType,local_corr

max_orbit_delay=300
max_clock_delay=30
# 因为Galileo的轨道和钟差可以一个历元解析出来，所以直接使用就行，不用这里个做复制，检查更新
class HASData:
    def __init__(self):
        self.init_empty()
    def init_empty(self):
        # 初始化所有属性为None或空列表/字典，保持结构一致但不包含实际数据
        self.cssrmode = None
        self.sat_n = []
        self.iodssr = None
        self.iodssr_c = {}
        self.lc = []
        self.nav_mode = {}
        self.subtype = None

    def update_value_from(self, source_object):
        # 从提供的对象复制属性，如果属性不存在则使用默认值
        self.cssrmode = getattr(source_object, 'cssrmode', None)
        self.sat_n = deepcopy(getattr(source_object, 'sat_n', []))
        self.iodssr = deepcopy(getattr(source_object, 'iodssr', None))
        self.iodssr_c = deepcopy(getattr(source_object, 'iodssr_c', []))
        self.nav_mode = deepcopy(getattr(source_object, 'nav_mode', []))
        self.subtype = deepcopy(getattr(source_object, 'subtype', None))
        self.mask_id=deepcopy(getattr(source_object, 'mask_id', None))
        self.mask_id_clk=deepcopy(getattr(source_object, 'mask_id_clk', None))
        self.sat_n_p=deepcopy(getattr(source_object, 'sat_n_p', None))
        # self.lc = deepcopy(getattr(source_object, 'lc', []))

        if len(self.lc)==0:
            self.lc.append(local_corr())
            inet = 0
            self.lc[inet].dclk = {}
            self.lc[inet].dorb = {}
            self.lc[inet].iode = {}
            self.lc[inet].iodc = {}
            self.lc[inet].iodc_c = {}
            self.lc[inet].cbias = {}
            self.lc[inet].t0={}
        # self.lc[0].iode=source_object.lc[0].iode
        self.lc[0].iode = deepcopy(getattr(source_object.lc[0], 'iode', {}))
        for j, sat in enumerate(source_object.sat_n):
            if source_object.iodssr >= 0 and source_object.iodssr_c[sCType.ORBIT] == source_object.iodssr:
                if sat not in source_object.sat_n:
                    continue
            if sat not in source_object.lc[0].iode.keys():
                continue
            if source_object.lc[0].dorb[sat] is None:
                continue
            if sat not in self.lc[0].t0:
                self.lc[0].t0[sat] = {}
            self.lc[0].iode[sat] = deepcopy(source_object.lc[0].iode[sat])
            self.lc[0].dorb[sat] = deepcopy(source_object.lc[0].dorb[sat])
            self.lc[0].t0[sat][sCType.ORBIT] = deepcopy(source_object.lc[0].t0[sat][sCType.ORBIT])
            if sat not in source_object.sat_n_p:
                print("missing clock corrections for sat="+str(sat))
                continue
            else:
                self.lc[0].dclk[sat] = deepcopy(source_object.lc[0].dclk[sat])
            self.lc[0].t0[sat][sCType.CLOCK] = deepcopy(source_object.lc[0].t0[sat][sCType.CLOCK])

    def deletePRN(self,sat):
        self.lc[0].iode[sat] = 0
        self.lc[0].dorb[sat] = []
        self.lc[0].t0[sat][sCType.ORBIT] = None

        self.lc[0].dclk[sat] = np.nan
        self.lc[0].t0[sat][sCType.CLOCK] = None

# Start epoch and number of epochs
#
start_date = datetime(2024, 5, 14)
process_days = 1
file_has_template = r'D:\work_lewen\source_code\git_lewen\NavDecoder\test_data\SEPT{}0.{}__SBF_GALRawCNAV.txt'
# nav_file_template = r'E:\GNSS_Data\products\eph\BRD400DLR_S_{}0000_01D_MN.rnx'
nav_file_template = r'D:\work_lewen\source_code\git_lewen\NavDecoder\test_data\BRDC00GOP_R_{}0000_01D_MN.rnx'
corr_dir_template = r'D:\work_lewen\source_code\git_lewen\NavDecoder\test_data\SEPT{}_HAS_new1'

for i in range(process_days):
    current_date = start_date + timedelta(days=i)
    ep = [current_date.year, current_date.month, current_date.day,
                  current_date.hour, current_date.minute, current_date.second]
    doy = current_date.timetuple().tm_yday
    year = current_date.year
    formatted_date = f"{year}{str(doy).zfill(3)}" 

    # 替换文件名中的年积日和年份
    file_has = file_has_template.format(str(doy).zfill(3),year-2000)
    nav_file = nav_file_template.format(formatted_date)

    # extend the navigation file to 3-days
    previous_date = current_date - timedelta(days=1)
    yyyy_doy0 = f"{year}{str(previous_date.timetuple().tm_yday).zfill(3)}" 
    nav_file0 = nav_file_template.format(yyyy_doy0)

    next_date = current_date + timedelta(days=1)
    yyyy_doy2 = f"{year}{str(next_date.timetuple().tm_yday).zfill(3)}" 
    nav_file2 = nav_file_template.format(yyyy_doy2)

    # generate the output file
    corr_dir = corr_dir_template.format(formatted_date)
    parent_dir = os.path.dirname(corr_dir)
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)
    print("=============Saving sp3/ssr/log to dir: "+corr_dir)
    file_sp3 = corr_dir + '.sp3'
    file_ssr = corr_dir + '.ssr'
    file_log = corr_dir + '.log'

    cs = cssr_has(file_log)
    cs.monlevel = 2
    time = epoch2time(ep)
    start_time = time
    week, tow = time2gpst(time)
    doy=time2doy(time)
    cs.week = week
    cs.tow0 = tow//86400*86400
    if not os.path.exists(file_has):
        continue

    # Read the raw HAS binary file according to the format of the Septentrio stardard
    dtype = [('tow', 'float64'), ('wn', 'int'),  ('prn', 'S3'), ('validity', 'S10'),('mask','int'),
             ('signal', 'S10'), ('num1', 'int'),('num2', 'int'), ('nav', 'S278')]
    v = np.genfromtxt(file_has, dtype=dtype, delimiter=',')
    v = v[v['validity']== b'Passed']
    # Eliminate whitespace
    for i, (nav, prn) in enumerate(zip(v['nav'], v['prn'])):
        v[i]['nav'] = b''.join(nav.split())
        v[i]['prn'] = int(prn[1:])

    # Read the Galileo-HAS Solomon matrix
    file_gm = r"B2b_HAS_decoder\Galileo-HAS-SIS-ICD_1.0_Annex_B_Reed_Solomon_Generator_Matrix.txt"
    gMat = np.genfromtxt(file_gm, dtype="u1", delimiter=",")

    # Read the navigation file from the successive 3-day
    rnx = rnxdec()
    nav = Nav()
    orb = peph()
    nav = rnx.decode_nav(nav_file0, nav)
    nav = rnx.decode_nav(nav_file, nav,True)
    nav = rnx.decode_nav(nav_file2, nav,True)
    nav_out = Nav()
    sp_out = peph()
    orb = peph()

    # 这两个变量是所有卫星通用的，如果某刻卫星中断后重新出现，基于这个时间是没办法判断的
    # 所以这个变量只能用来判断数据是否历元更新，因为历元是按照时间来的
    record_orbit_update_time=None
    record_clock_update_time=None
    orbit_data={}
    clock_data={}
    HASData0=HASData()
    delay=0

    mid_ = -1
    ms_ = -1
    icnt = 0
    rec = []
    mid_decoded = []
    has_pages = np.zeros((255, 53), dtype=int)
    valid_tows = np.unique(v['tow'])
    current_time=start_time
    for hh in tqdm(range(len(valid_tows))):
        decode_page=False
        tow = valid_tows[hh]
        cs.tow0 = tow // 3600 * 3600
        vi = v[v['tow'] == tow]
        for vn in vi:
            buff = unhexlify(vn['nav'])
            i = 14
            if bs.unpack_from('u24', buff, i)[0] == 0xaf3bc3:
                continue
            hass, res = bs.unpack_from('u2u2', buff, i)
            i += 4
            if hass >= 2:  # 0:test,1:operational,2:res,3:dnu
                continue
            mt, mid, ms, pid = bs.unpack_from('u2u5u5u8', buff, i)

            cs.msgtype = mt
            ms += 1
            i += 20

            if mid_ == -1 and mid not in mid_decoded:
                mid_ = mid
                ms_ = ms
            if mid == mid_ and pid-1 not in rec:
                page = bs.unpack_from('u8'*53, buff, i)
                rec += [pid-1]
                has_pages[pid-1, :] = page

            # print(f"{mt} {mid} {ms} {pid}")

        if len(rec) >= ms_:
            if cs.monlevel >= 2:
                print("data collected mid={:2d} ms={:2d} tow={:.0f}"
                      .format(mid_, ms_, tow))
            HASmsg = cs.decode_has_page(rec, has_pages, gMat, ms_)
            cs.decode_cssr(HASmsg)
            decode_page=True
            rec = []
            mid_decoded += [mid_]
            mid_ = -1
            if len(mid_decoded) > 10:
                mid_decoded = mid_decoded[1:]
        else:
            icnt += 1
            if icnt > 2*ms_ and mid_ != -1:
                icnt = 0
                if cs.monlevel >= 2:
                    print(f"reset mid={mid_} ms={ms_} tow={tow}")
                rec = []
                mid_ = -1
        intervals=5
        if decode_page:
            update_Orbssr=False
            update_Clkssr=False
            newssr_time=None
            lastssr_time=None
            if cs.subtype == sCSSR.ORBIT or cs.subtype==sCSSR.CBIAS:
                if record_orbit_update_time is None:
                    record_orbit_update_time=cs.time
                time_orbit_sat=cs.time
                if abs(timediff(time_orbit_sat, record_orbit_update_time))>1: #comes new orbit
                    update_Orbssr=True
                    newssr_time=time_orbit_sat
                    lastssr_time=record_orbit_update_time
            if cs.subtype == sCSSR.CLOCK:
                if record_clock_update_time is None:
                    record_clock_update_time = cs.time
                time_clock_sat=cs.time
                if timediff(time_clock_sat, record_clock_update_time)>=1: #comes new clock
                    update_Clkssr=True
                    newssr_time=time_clock_sat
                    lastssr_time=record_clock_update_time
            if update_Clkssr or update_Orbssr:
                str_obs1 = time2str(time_clock_sat)
                str_obs2 = time2str(record_clock_update_time)
                time_debug = epoch2time([2024, 3, 16, 0, 0, 45])
                # 根据current_time查找最新，可用的产品，实时的产品时间应远于目前的
                str_obs=time2str(current_time)
                if abs(timediff(current_time, time_debug))<1:
                    print(time2str(time_debug))
                while timediff(current_time,newssr_time)<0:
                    time_corr = timeadd(current_time, -delay)
                    debug_obs=time2str(current_time)
                    if timediff(time_corr, lastssr_time)<=0:
                        current_time=timeadd(time_corr, intervals)
                        continue
                    if update_Clkssr and timediff(time_corr, lastssr_time)>max_clock_delay:
                        cs.log_msg(">>>>ERROR: large clock difference[obst-clkt] : " + time2str(time_corr) + " " + time2str(lastssr_time))
                        current_time=timeadd(time_corr, intervals)
                        continue
                    if update_Orbssr and timediff(time_corr, record_orbit_update_time)<0 or timediff(time_corr, record_orbit_update_time)>max_orbit_delay:
                        cs.log_msg(">>>>ERROR: large orbit difference [obst-orbt]: " + time2str(time_corr) + " " + time2str(record_orbit_update_time))
                        current_time=timeadd(time_corr, intervals)
                        continue
                    # cs.encode_SP3(HASData0,orb, nav, current_time, record_clock_update_time,sp_out, nav_out, file_ssr)
                    cs.encode_SP3(HASData0,orb, nav, current_time, sp_out, nav_out, file_ssr)
                    current_time=timeadd(time_corr, intervals)
                if update_Clkssr:
                    record_clock_update_time=time_clock_sat
                if update_Orbssr:
                    record_orbit_update_time = time_orbit_sat
                if cs.mask_id ==cs.mask_id_clk:
                    HASData0.update_value_from(cs)
    sp_out.write_sp3(file_sp3, nav_out)