"""
 Decoding the PPP-B2b corrections from Septentrio receiver
"""
from binascii import unhexlify
from copy import deepcopy
import numpy as np
import os
from B2b_HAS_decoder.gnss import *
from B2b_HAS_decoder.peph import peph
from B2b_HAS_decoder.cssr_bds_sept import cssr_bds
from B2b_HAS_decoder.rinex import rnxdec
from B2b_HAS_decoder.sdr_ldpc_test import *
from B2b_HAS_decoder.cssrlib import sCSSR,sCType,local_corr
from datetime import datetime, timedelta

class B2BData:
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
        self.iodssr = getattr(source_object, 'iodssr', None)
        self.iodssr_c = deepcopy(getattr(source_object, 'iodssr_c', []))
        self.nav_mode = deepcopy(getattr(source_object, 'nav_mode', []))
        self.subtype = getattr(source_object, 'subtype', None)
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
            # 钟差和轨道IOD匹配，使用新的钟差
            if source_object.lc[0].iodc[sat] == source_object.lc[0].iodc_c[sat]:
                self.lc[0].iode[sat] = deepcopy(source_object.lc[0].iode[sat])
                self.lc[0].dorb[sat] = deepcopy(source_object.lc[0].dorb[sat])
                self.lc[0].iodc[sat] = deepcopy(source_object.lc[0].iodc[sat])
                self.lc[0].t0[sat][sCType.ORBIT] = deepcopy(source_object.lc[0].t0[sat][sCType.ORBIT])

                self.lc[0].dclk[sat] = deepcopy(source_object.lc[0].dclk[sat])
                self.lc[0].iodc_c[sat] = deepcopy(source_object.lc[0].iodc_c[sat])
                self.lc[0].t0[sat][sCType.CLOCK] = deepcopy(source_object.lc[0].t0[sat][sCType.CLOCK])
            # 不更新钟差的改正数和IOD
            else:
                self.lc[0].iode[sat] = deepcopy(source_object.lc[0].iode[sat])
                self.lc[0].dorb[sat] = deepcopy(source_object.lc[0].dorb[sat])
                self.lc[0].iodc[sat] = deepcopy(source_object.lc[0].iodc[sat])
                self.lc[0].t0[sat][sCType.ORBIT] = deepcopy(source_object.lc[0].t0[sat][sCType.ORBIT])

    def deletePRN(self,sat):
        self.lc[0].iode[sat] = 0
        self.lc[0].dorb[sat] = []
        self.lc[0].iodc[sat] = 0
        self.lc[0].t0[sat][sCType.ORBIT] = None

        self.lc[0].dclk[sat] = np.nan
        self.lc[0].iodc_c[sat] = 0
        self.lc[0].t0[sat][sCType.CLOCK] = None

# Main setting for the processing
start_date = datetime(2024, 5, 14)
process_days = 1

max_orbit_delay=300
max_clock_delay=30
file_bds_template = r'D:\work_lewen\source_code\git_lewen\NavDecoder\test_data\SEPT{}0.{}__SBF_BDSRawB2b.txt'
nav_file_template = r'D:\work_lewen\source_code\git_lewen\NavDecoder\test_data\BRD400DLR_S_{}0000_01D_MN.rnx'
corr_dir_template = r'D:\work_lewen\source_code\git_lewen\NavDecoder\test_data\SEPT{}_B2B_new'
for i in range(process_days):
    current_date = start_date + timedelta(days=i)
    ep = [current_date.year, current_date.month, current_date.day,
                  current_date.hour, current_date.minute, current_date.second]
    doy = current_date.timetuple().tm_yday
    year = current_date.year
    formatted_date = f"{year}{str(doy).zfill(3)}" 

    file_bds = file_bds_template.format(str(doy).zfill(3),year-2000)
    nav_file = nav_file_template.format(formatted_date)

    if not os.path.exists(file_bds):
        print("File not found: "+file_bds)
        continue

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
    print("=============Saving sp3/ssr/log to dir: " + corr_dir)
    file_sp3 = corr_dir + '.sp3'
    file_ssr = corr_dir + '.ssr'
    file_log = corr_dir + '.log'
    cs = cssr_bds(file_log)
    cs.monlevel = 2
    prn_ref = 59  # satellite PRN to receive BDS PPP collection

    time = epoch2time(ep)
    current_time = time
    week, tow = time2gpst(time)
    doy=time2doy(time)
    cs.week = week
    cs.tow0 = tow//86400*86400
    # Read the PPP-B2b binary file
    dtype = [('tow', 'float64'), ('wn', 'int'),  ('prn', 'S3'), ('validity', 'S10'),
             ('signal', 'S10'), ('num2', 'int'), ('nav', 'S278')]
    v = np.genfromtxt(file_bds, dtype=dtype, delimiter=',')
    v = v[v['validity']== b'Passed']
    v = v[v['prn'] == ('C' + str(prn_ref)).encode()]
    # Eliminate whitespace
    for i, (nav, prn) in enumerate(zip(v['nav'], v['prn'])):
        v[i]['nav'] = b''.join(nav.split())
        v[i]['prn'] = int(prn[1:])

    # Read the navigation file
    rnx = rnxdec()
    nav = Nav()
    orb = peph()
    nav = rnx.decode_nav(nav_file, nav)
    nav = rnx.decode_nav(nav_file, nav,True)
    nav = rnx.decode_nav(nav_file2, nav,True)
    nav_out = Nav()
    sp_out = peph()

    record_orbit_update_time=None
    record_clock_update_time=None
    orbit_data={}
    clock_data={}
    B2BData0=B2BData()
    delay=0
    for row in v:
        buff = decode_LDPC(row)
        mt=cs.decode_cssr(buff, 0)
        intervals=5
        if (cs.lc[0].cstat & 0xf) == 0xf:
            if cs.subtype == sCSSR.CLOCK:
                if record_clock_update_time is None:
                    record_clock_update_time = cs.time
                time_clock_sat=cs.time
                str_obs1 = time2str(time_clock_sat)
                str_obs2 = time2str(record_clock_update_time)
                time_test = epoch2time([2023, 12, 3, 0, 1, 55])
                if timediff(time_clock_sat, record_clock_update_time)>=1: #到这里，record_clock_update_time这一时刻的钟差已经全部解析完毕，可用
                    '''生成处理GNSS的时间间隔，为了模拟实时要求，保证观测时间要早于可用的轨道和钟差时间'''
                    # 根据current_time查找最新，可用的产品，实时的产品时间应远于目前的
                    str_obs=time2str(current_time)
                    if abs(timediff(current_time,time_test))<1:
                        print(time2str(time_test))
                    while timediff(current_time,time_clock_sat)<0:
                        time_corr = timeadd(current_time, -delay)
                        debug_obs=time2str(current_time)
                        if timediff(time_corr, record_clock_update_time)<0:
                            current_time=timeadd(time_corr, intervals)
                            continue
                        if timediff(time_corr, record_clock_update_time)>max_clock_delay:
                            cs.log_msg(">>>>ERROR: large clock difference[obst-clkt] : " + time2str(time_corr) + " " + time2str(record_clock_update_time))
                            current_time=timeadd(time_corr, intervals)
                            continue
                        if timediff(time_corr, record_orbit_update_time)<0 or timediff(time_corr, record_orbit_update_time)>max_orbit_delay:
                            cs.log_msg(">>>>ERROR: large orbit difference [obst-orbt]: " + time2str(time_corr) + " " + time2str(record_orbit_update_time))
                            current_time=timeadd(time_corr, intervals)
                            continue
                        cs.encode_SP3(B2BData0, orb, nav, current_time, record_clock_update_time,sp_out, nav_out, file_ssr)
                        current_time=timeadd(time_corr, intervals)
                    record_clock_update_time=time_clock_sat
                else:
                    B2BData0.update_value_from(cs)

            if cs.subtype == sCSSR.ORBIT:
                if record_orbit_update_time is None:
                    record_orbit_update_time=cs.time
                time_orbit_sat=cs.time
                if abs(timediff(time_orbit_sat, record_orbit_update_time))>1:
                    '''这里轨道更新了，可能也就意味着进行广播星历匹配的IOD更新了，且之后解码得到的都是新的历元的钟差信息，所以这里也重置B2BData0保存的数据信息'''
                    record_orbit_update_time=time_orbit_sat
                    B2BData0.update_value_from(cs)

    sp_out.write_sp3(file_sp3, nav_out)