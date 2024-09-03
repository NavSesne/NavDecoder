"""
BDS PPP correction data decoder

[1] BeiDou Navigation Satellite System Signal In Space Interface Control Document
Precise Point Positioning Service Signal PPP-B2b (Version 1.0), 2020
"""

import numpy as np
import bitstruct as bs
from B2b_HAS_decoder.cssrlib import cssr, sCSSR, sCSSRTYPE, sGNSS, prn2sat, sCType
from B2b_HAS_decoder.gnss import bdt2time, bdt2gpst, uGNSS, uSIG, uTYP, rSigRnx,sat2prn,time2str,sat2id,time2epoch,timediff,vnorm,rCST
from B2b_HAS_decoder.ephemeris import eph2pos, findeph,eph2rel
from B2b_HAS_decoder.peph import peph,peph_t


max_orbit_delay=300
max_clock_delay=30

class cssr_bds(cssr):
    def __init__(self, foutname=None):
        super().__init__(foutname)
        self.MAXNET = 1
        self.cssrmode = sCSSRTYPE.BDS_PPP
        self.nsig_max = 8
        self.iodp = -1
        self.iodp_p = -1
        self.cb_blen = 12
        self.cb_scl = 0.017

    def ssig2rsig(self, sys: sGNSS, utyp: uTYP, ssig):
        gps_tbl = {
            0: uSIG.L1C,
            1: uSIG.L1P,
            4: uSIG.L1L,
            5: uSIG.L1X,
            7: uSIG.L2L,
            8: uSIG.L2X,
            11: uSIG.L5I,
            12: uSIG.L5Q,
            13: uSIG.L5X,
        }
        glo_tbl = {
            0: uSIG.L1C,
            1: uSIG.L1P,
            2: uSIG.L2C,
        }

        gal_tbl = {
            1: uSIG.L1B,
            2: uSIG.L1C,
            4: uSIG.L5Q,
            5: uSIG.L5I,
            7: uSIG.L7I,
            8: uSIG.L7Q,
            11: uSIG.L6C,
        }

        bds_tbl = {
            0: uSIG.L2I,
            1: uSIG.L1D,
            2: uSIG.L1P,
            4: uSIG.L5D,
            5: uSIG.L5P,
            7: uSIG.L7D,
            8: uSIG.L7P,
            12: uSIG.L6I,
        }

        usig_tbl_ = {
            uGNSS.GPS: gps_tbl,
            uGNSS.GLO: glo_tbl,
            uGNSS.GAL: gal_tbl,
            uGNSS.BDS: bds_tbl,
        }

        usig_tbl = usig_tbl_[sys]
        return rSigRnx(sys, utyp, usig_tbl[ssig])

    def sval(self, u, n, scl):
        """ calculate signed value based on n-bit int, lsb """
        invalid = -2**(n-1)
        dnu = -(2**(n-1)-1)  # this value seems to be invalid
        y = np.nan if u == invalid or u == dnu else u*scl
        return y

    def slot2prn(self, slot):
        prn = 0
        sys = uGNSS.NONE
        if slot >= 1 and slot <= 63:
            prn = slot
            sys = uGNSS.BDS
        elif slot <= 100:
            prn = slot-63
            sys = uGNSS.GPS
        elif slot <= 137:
            prn = slot-100
            sys = uGNSS.GAL
        elif slot <= 174:
            prn = slot-137
            sys = uGNSS.GLO
        return sys, prn

    def decode_head(self, msg, i, st=-1):
        self.tod, _, iodssr = bs.unpack_from('u17u4u2', msg, i)
        i += 23

        if st == sCSSR.MASK:
            self.iodssr = iodssr

        if self.tow0 >= 0:
            self.tow = self.tow0+self.tod
            if self.week >= 0:
                self.time = bdt2gpst(bdt2time(self.week, self.tow))

        head = {'uint': 0, 'mi': 0, 'iodssr': iodssr}
        return head, i

    def add_gnss(self, mask, blen, gnss):
        prn, nsat = self.decode_mask(mask, blen)
        self.nsat_g[gnss] = nsat
        self.nsat_n += nsat
        if nsat > 0:
            self.ngnss += 1
        sys = self.gnss2sys(gnss)
        for k in range(0, nsat):
            # 这里的sat是按照GPS,GAL,QZSS,GLO,BDS的顺序叠加的序号
            sat = prn2sat(sys, prn[k])
            # 对于每一颗卫星，都会保存期导航系统，卫星编号和
            self.sys_n.append(sys)
            self.sat_n.append(sat)
            self.gnss_n.append(gnss)

    def decode_cssr_mask(self, msg, i):
        """decode MT1 Mask message """
        head, i = self.decode_head(msg, i, sCSSR.MASK)
        # 这个数据很重要，只有当轨道钟差中的内容和这一致的时候才能进行使用
        self.iodp = bs.unpack_from('u4', msg, i)[0]
        i += 4

        mask_bds, mask_gps, mask_gal, mask_glo = \
            bs.unpack_from('u63u37u37u37', msg, i)
        i += 174

        if self.iodp != self.iodp_p: #iodp_p表示之前历元的iodp
            # 这里是对之前变量的值进行初始化或者更新
            self.sat_n_p = self.sat_n
            self.iodssr_p = self.iodssr
            self.sig_n_p = self.sig_n
            self.iodp_p = self.iodp

            self.nsat_n = 0 #解码得到的播发改正数的卫星总数量，包含所有的导航系统，比如为59
            self.nsig_n = []
            self.sys_n = []#对应上面59颗卫星中每颗卫星的导航系统
            self.gnss_n = []#对应上面59颗卫星中每颗卫星的导航系统，和sys_n基本一样，两者通过gnss2sys函数转换
            self.sat_n = []#对应上面59颗卫星中每颗卫星的编号（从GPS的1号星开始的累计编号）
            self.nsig_total = 0
            self.sig_n = []
            self.nm_idx = np.zeros(self.SYSMAX, dtype=int)
            self.ngnss = 0 #可用的导航系统数（如GPS+BDS，该值为2）
            # self.gnss_idx = np.zeros(self.ngnss, dtype=int)
            self.nsat_g = np.zeros(self.SYSMAX, dtype=int)

            self.add_gnss(mask_bds, 63, sGNSS.BDS)
            self.add_gnss(mask_gps, 37, sGNSS.GPS)
            self.add_gnss(mask_gal, 37, sGNSS.GAL)
            self.add_gnss(mask_glo, 37, sGNSS.GLO)

            # 这里是针对每颗卫星进行存储器轨道、钟差的改正数变量进行初始化
            inet = 0
            self.lc[inet].dclk = {}
            self.lc[inet].dorb = {}
            self.lc[inet].iode = {}
            self.lc[inet].iodc = {}
            self.lc[inet].iodc_c = {}
            self.lc[inet].cbias = {}
            self.nsig_n = np.ones(self.nsat_n, dtype=int)*self.nsig_max
            self.sig_n = {}
            self.ura = {}

            # fallback for inconsistent clock update
            self.lc[inet].dclk_p = {}
            self.lc[inet].iodc_c_p = {}
            msg="Change of IODP in decode_cssr_mask"
            self.log_msg(msg)
        # 这里存储的是对应于参数解码中的IOD SSR (SSR版本号）
        # IOD SSR 变化表明数据生成配置的变化。各类型数据首先需要保证各自 IOD SSR 相同，才可匹配使用
        self.iodssr = head['iodssr']

        self.lc[0].cstat |= (1 << sCType.MASK)
        # self.lc[0].t0[0][sCType.MASK] = self.time
        self.log_msg("decode_ssr_mask: iodp= "+str(self.iodp))
        return i

    def decode_cssr_orb_sat(self, msg, i, inet, sat_n):
        slot, iodn, iodc = bs.unpack_from('u9u10u3', msg, i)
        i += 22

        dx, dy, dz = bs.unpack_from('s{:d}s{:d}s{:d}'.format(
            self.dorb_blen[0], self.dorb_blen[1], self.dorb_blen[2]), msg, i)
        i += self.dorb_blen[0] + self.dorb_blen[1] + self.dorb_blen[2]

        ucls, uval = bs.unpack_from('u3u3', msg, i)
        i += 6

        sys, prn = self.slot2prn(slot)
        if sys == uGNSS.NONE:
            return i
        sat = prn2sat(sys, prn)
        if sat not in sat_n:
            return i

        if (sys == uGNSS.GPS) or (sys == uGNSS.BDS):
            iodn = iodn & 0xff  # IODC -> IODE

        dorb = np.zeros(3)
        dorb[0] = self.sval(dx, self.dorb_blen[0], self.dorb_scl[0])
        dorb[1] = self.sval(dy, self.dorb_blen[1], self.dorb_scl[1])
        dorb[2] = self.sval(dz, self.dorb_blen[2], self.dorb_scl[2])

        self.lc[inet].iode[sat] = iodn
        self.lc[inet].iodc[sat] = iodc
        self.lc[inet].dorb[sat] = dorb
        self.ura[sat] = self.quality_idx(ucls, uval)

        if sat not in self.lc[inet].t0:
            self.lc[inet].t0[sat] = {}

        self.lc[inet].t0[sat][sCType.ORBIT] = self.time

        str_cssr_orb_sat = 'cssr_orb_sat: %s %s %d %.3f %.3f %.3f ' %(sat2id(sat), time2str(self.time), iodc, dorb[0], dorb[1], dorb[2])
        self.log_msg(str_cssr_orb_sat)

        return i

    def decode_cssr_orb(self, msg, i, inet=0):
        """decode MT2 orbit + URA message """
        head, i = self.decode_head(msg, i)
        sat_n = np.array(self.sat_n)

        if self.iodssr != head['iodssr']:
            return -1

        # 解码得到的变量存储在下面的数组里
        # self.lc[inet].iode[sat] = iodn
        # self.lc[inet].iodc[sat] = iodc
        # self.lc[inet].dorb[sat] = dorb
        # self.iodssr_c[sCType.ORBIT]
        # self.lc[inet].t0[sat][sCType.ORBIT] = self.time
        for k in range(6):
            i = self.decode_cssr_orb_sat(msg, i, inet, sat_n)

        self.iodssr_c[sCType.ORBIT] = head['iodssr']
        self.lc[inet].cstat |= (1 << sCType.ORBIT)

        i += 19  # reserve
        return i

    def decode_cssr_cbias(self, msg, i, inet=0):
        """decode MT3 Code Bias Correction message """
        head, i = self.decode_head(msg, i)
        nsat = self.nsat_n
        self.flg_net = False
        if self.iodssr != head['iodssr']:
            return -1

        sat_n = np.array(self.sat_n)
        nsat = bs.unpack_from('u5', msg, i)[0]
        i += 5

        for k in range(nsat):
            slot, nsig = bs.unpack_from('u9u4', msg, i)
            i += 13
            sys, prn = self.slot2prn(slot)
            sat = prn2sat(sys, prn)
            if sat not in sat_n:
                continue
            if sat not in self.lc[inet].t0:
                self.lc[inet].t0[sat] = {}
            self.lc[inet].t0[sat][sCType.CBIAS] = self.time

            self.lc[inet].cbias[sat] = {}
            for j in range(0, nsig):
                sig, cb = bs.unpack_from(
                    'u4s{:d}'.format(self.cb_blen), msg, i)
                i += 4+self.cb_blen
                sig_ = self.ssig2rsig(sys, uTYP.C, sig)
                self.lc[inet].cbias[sat][sig_] = self.sval(
                    cb, self.cb_blen, self.cb_scl)

        self.iodssr_c[sCType.CBIAS] = head['iodssr']
        self.lc[inet].cstat |= (1 << sCType.CBIAS)

        return i

    def decode_cssr_clk_sat(self, msg, i, inet, sat):
        """ decode clock correction for satellite """
        iodc, dclk = bs.unpack_from('u3s{:d}'.format(self.dclk_blen), msg, i)
        i += 3+self.dclk_blen

        if sat in self.lc[inet].iodc_c.keys() and \
                iodc != self.lc[inet].iodc_c[sat]:
            self.lc[inet].iodc_c_p[sat] = self.lc[inet].iodc_c[sat]
            self.lc[inet].dclk_p[sat] = self.lc[inet].dclk[sat]
        self.lc[inet].iodc_c[sat] = iodc
        # note: the sign of the clock correction reversed
        self.lc[inet].dclk[sat] = - \
            self.sval(dclk, self.dclk_blen, self.dclk_scl)
        if ~np.isnan(self.lc[inet].dclk[sat]):
            str_cssr_clk_sat = 'cssr_clk_sat: %s %s %d %.3f ' % (sat2id(sat), time2str(self.time), iodc, self.lc[inet].dclk[sat])
            self.log_msg(str_cssr_clk_sat)
        return i

    def decode_cssr_clk(self, msg, i, inet=0):
        """decode MT4 Clock Correction message """
        head, i = self.decode_head(msg, i)
        if self.iodssr != head['iodssr']:
            return -1

        iodp, st1 = bs.unpack_from('u4u5', msg, i)
        i += 9
        # 这里的self.iodp就是mask解码得到的，如果不匹配，就跳过对应内容的解码
        if iodp != self.iodp:
            i += 23*18+10
            return i
        # 如果iodp匹配，则解码23颗卫星钟差的钟差参数,得到：
        # 卫星钟差的IODC:   self.lc[inet].iodc_c[sat] = iodc 这个每颗卫星不一样
        # 卫星钟差改正：     self.lc[inet].dclk[sat] =
        # 这颗卫星的时间，   self.lc[inet].t0[sat][sCType.CLOCK] = self.time 这里分别保存了每一种参数的时间
        # 钟差的iod_ssr:   self.iodssr_c[sCType.CLOCK] = head['iodssr']
        for k in range(23):
            idx = st1*23+k
            if idx < self.nsat_n:
                sat = self.sat_n[idx]
                i = self.decode_cssr_clk_sat(msg, i, inet, sat)
                if sat not in self.lc[inet].t0:
                    self.lc[inet].t0[sat] = {}
                self.lc[inet].t0[sat][sCType.CLOCK] = self.time

        self.iodssr_c[sCType.CLOCK] = head['iodssr']
        self.lc[inet].cstat |= (1 << sCType.CLOCK)

        i += 10
        return i

    def decode_cssr_ura(self, msg, i):
        """decode MT5 URA message """
        head, i = self.decode_head(msg, i)
        if self.iodssr != head['iodssr']:
            return -1

        iodp, st2 = bs.unpack_from('u4u3', msg, i)
        i += 7

        for k in range(70):
            v = bs.unpack_from_dict('u3u3', ['class', 'val'], msg, i)
            i += 6
            idx = st2*70+k
            if idx < self.nsat_n:
                sat = self.sat_n[idx]
                self.ura[sat] = self.quality_idx(v['class'], v['val'])

        self.lc[0].cstat |= (1 << sCType.URA)
        # self.lc[0].t0[0][sCType.URA] = self.time
        return i

    def decode_cssr_comb1(self, msg, i, inet=0):
        """decode MT6 combined message #1 """
        numc, numo = bs.unpack_from('u5u3', msg, i)
        sat_n = np.array(self.sat_n)

        if numc > 0:
            tod, _, iodssr, iodp, slot_s = \
                bs.unpack_from('u17u4u2u4u9', msg, i)
            i += 36
            for k in range(numc):
                idx = slot_s+k
                sat = self.sat_n[idx]
                i = self.decode_cssr_clk_sat(msg, i, inet, sat)

        if numo > 0:
            tod, _, iodssr = bs.unpack_from('u17u4u2', msg, i)
            i += 23
            for k in range(numo):
                i = self.decode_cssr_orb_sat(msg, i, inet, sat_n)

        return i

    def decode_cssr_comb2(self, msg, i, inet=0):
        """decode MT7 combined message #2 """
        numc, numo = bs.unpack_from('u5u3', msg, i)
        sat_n = np.array(self.sat_n)

        if numc > 0:
            tod, _, iodssr = bs.unpack_from('u17u4u2', msg, i)
            i += 23
            for k in range(numc):
                slot = bs.unpack_from('u9', msg, i)
                i += 9
                sys, prn = self.slot2prn(slot)
                sat = prn2sat(sys, prn)
                i = self.decode_cssr_clk_sat(msg, i, inet, sat)

        if numo > 0:
            tod, _, iodssr = bs.unpack_from('u17u4u2', msg, i)
            i += 23
            for k in range(numo):
                i = self.decode_cssr_orb_sat(msg, i, inet, sat_n)

        return i

    def decode_cssr(self, msg, i=0):
        mt = bs.unpack_from('u6', msg, i)[0]
        i += 6

        if mt == 1:
            self.subtype = sCSSR.MASK
            i = self.decode_cssr_mask(msg, i)
        elif mt == 2:
            self.subtype = sCSSR.ORBIT
            i = self.decode_cssr_orb(msg, i)
            # if self.monlevel > 0 and self.fh is not None:
            #     self.out_log()
        elif mt == 3:
            self.subtype = sCSSR.CBIAS
            i = self.decode_cssr_cbias(msg, i)
            if self.monlevel > 0 and self.fh is not None:
                self.out_log()
        elif mt == 4:
            self.subtype = sCSSR.CLOCK
            i = self.decode_cssr_clk(msg, i)
            # if self.monlevel > 0 and self.fh is not None:
            #     self.out_log()
        elif mt == 5:
            self.subtype = sCSSR.URA
            i = self.decode_cssr_ura(msg, i)
        elif mt == 6:
            i = self.decode_cssr_comb1(msg, i)
        elif mt == 7:
            i = self.decode_cssr_comb2(msg, i)

        # msg=" mt={:2d} tow={:6.1f} \n".format(mt, self.tow)
        # self.log_msg(msg)

        return mt

    def log_msg(self,msg):
        if self.monlevel > 0 and self.fh is not None:
            self.fh.write(msg+"\n")
    def encode_BNC_corr(self,row,run_orbit_update_time,run_clock_update_time):
        file_ssr=r"D:\SynologyDrive\202311\log_ssr.txt"
        orbit_str=[]
        clock_str=[]
        time_epoch_orbit=None
        time_epoch_clock=None
        """首先进行IOD匹配和卫星匹配检查"""
        for k in range(uGNSS.MAXSAT):
            sat = k+1
            sys, prn = sat2prn(sat)
            sat_id= sat2id(sat)
            if not (sys == uGNSS.GPS or sys==uGNSS.BDS):
                continue
            # cs.iodssr: mask里面解码得到的
            #iodssr_c[sCType.ORBIT] 是轨道部分解码得到的IOD
            if self.iodssr >= 0 and self.iodssr_c[sCType.ORBIT] == self.iodssr:
                # 判断解算的卫星是否在B2b解码的改正数里面
                if sat not in self.sat_n:
                    continue
            # 本历元的不一样，则轨道是否与判断之前历元的IOD一致
            elif self.iodssr_p >= 0 and self.iodssr_c[sCType.ORBIT] == self.iodssr_p:
                if sat not in self.sat_n_p:
                    continue
            else:
                continue
            # 判断这颗卫星是否有轨道改正
            if sat not in self.lc[0].iode.keys():
                continue

            iode = self.lc[0].iode[sat]
            dorb = self.lc[0].dorb[sat]  # radial,along-track,cross-track

            if self.cssrmode == sCSSRTYPE.BDS_PPP:  # consistency check for IOD corr
                # cs.lc[0].iodc[sat] 是轨道里面对应的IODC
                # cs.lc[0].iodc_c[sat] 是钟差里保存的IODC
                if self.lc[0].iodc[sat] == self.lc[0].iodc_c[sat]:
                    dclk = self.lc[0].dclk[sat]
                else:
                    if self.lc[0].iodc[sat] == self.lc[0].iodc_c_p[sat]:
                        dclk = self.lc[0].dclk_p[sat]
                    else:
                        continue

            if np.isnan(dclk) or np.isnan(dorb@dorb):
                continue
            if sat in self.lc[0].t0 and sCType.ORBIT in self.lc[0].t0[sat]:
                if run_orbit_update_time is None:
                    run_orbit_update_time= self.lc[0].t0[sat][sCType.ORBIT]
            if sat in self.lc[0].t0 and sCType.CLOCK in self.lc[0].t0[sat]:
                if run_clock_update_time is None:
                    run_clock_update_time = self.lc[0].t0[sat][sCType.CLOCK]
            time_orbit_sat=self.lc[0].t0[sat][sCType.ORBIT]
            time_clock_sat=self.lc[0].t0[sat][sCType.CLOCK]
            if timediff(time_orbit_sat,run_orbit_update_time)>1:
                continue
            else:
                # 符合BNC格式的输出
                # str_orbit=('{}  {:4d} {:11.4f} {:11.4f} {:11.4f} \n'.format(sat_id, self.lc[0].iodc[sat], dorb[0], dorb[1], dorb[2]))
                '''以下是测试格式的输出'''
                str_orbit=('{}  {:4d} {:4d} {:11.4f} {:11.4f} {:11.4f} \n'.format(sat_id, self.iodssr, self.lc[0].iodc[sat], dorb[0], dorb[1], dorb[2]))
                orbit_str.append(str_orbit)
            if timediff(time_clock_sat,run_clock_update_time)>1:
                continue
            else:
                # 符合BNC格式的输出
                # str_clock=('{}  {:4d} {:11.4f} {:11.4f} {:11.4f} \n'.format(sat_id, self.lc[0].iodc_c[sat],dclk,0.0,0.0))
                '''以下是测试格式的输出'''
                # self.iodssr: MASK解码得到的IOD
                # self.lc[inet].iodc_c[sat] = iodc: 卫星钟差的IODC,这个每颗卫星不一样
                # self.iodssr_c[sCType.CLOCK]:      帧头里的SSR，同一帧数据所有卫星一致
                str_clock=('{}  {:4d} {:4d} {:4d} {:11.4f} {:11.4f} {:11.4f} \n'.format(sat_id, self.iodssr, self.iodssr_c[sCType.CLOCK],self.lc[0].iodc_c[sat],dclk,0.0,0.0))
                clock_str.append(str_clock)
        nsat_orbit=len(orbit_str)
        e = time2epoch(time_epoch_orbit)
        orbit_time="{:04d} {:02d} {:02d} {:02d} {:02d} {:02d}".format(e[0], e[1], e[2], e[3], e[4], int(e[5]))
        e = time2epoch(time_epoch_clock)
        clock_time="{:04d} {:02d} {:02d} {:02d} {:02d} {:02d}".format(e[0], e[1], e[2], e[3], e[4], int(e[5]))
        corr_str = np.array2string(row, separator=', ')
        with open(file_ssr,'a') as fp:
            fp.write(corr_str+"\n")
            fp.write('> ORBIT {} {}  {:4d} {} \n'.format(orbit_time, 2, nsat_orbit,"B2B"))
            for ln in orbit_str:
                fp.write(ln)
            fp.write('> CLOCK {} {}  {:4d} {} \n'.format(clock_time, 2, nsat_orbit,"B2B"))
            for ln in clock_str:
                fp.write(ln)

    def encode_BNC_corr_new(self,row,orbit_data,run_orbit_update_time,clock_data,run_clock_update_time):
        file_ssr=r"D:\SynologyDrive\202311\log_ssr.txt"
        if self.cssrmode != sCSSRTYPE.BDS_PPP:  # consistency check for IOD corr
            return

        """首先进行IOD匹配和卫星匹配检查"""
        for k in range(uGNSS.MAXSAT):
            sat = k+1
            sys, prn = sat2prn(sat)
            sat_id= sat2id(sat)
            if not (sys == uGNSS.GPS or sys==uGNSS.BDS):
                continue
            # cs.iodssr: mask里面解码得到的
            #iodssr_c[sCType.ORBIT] 是轨道部分解码得到的IOD
            if self.iodssr >= 0 and self.iodssr_c[sCType.ORBIT] == self.iodssr:
                # 判断解算的卫星是否在B2b解码的改正数里面
                if sat not in self.sat_n:
                    continue
            # 判断这颗卫星是否有轨道改正
            if sat not in self.lc[0].iode.keys():
                continue
            if self.subtype == sCSSR.CLOCK:
                iode = self.lc[0].iode[sat]
                if self.lc[0].iodc[sat] == self.lc[0].iodc_c[sat]:
                    dclk = self.lc[0].dclk[sat]
                else:
                    print("Different IOD definition: " + str(self.lc[0].iodc[sat]) + " " + str(self.lc[0].iodc_c[sat]))
                    continue
                if np.isnan(dclk):
                    continue
                if run_clock_update_time is None:
                    run_clock_update_time = self.lc[0].t0[sat][sCType.CLOCK]
                time_clock_sat=self.lc[0].t0[sat][sCType.CLOCK]
                if timediff(time_clock_sat,run_clock_update_time)>1:
                    # epoch update, output the products
                    nsat_clock=len(clock_data)
                    e = time2epoch(run_clock_update_time)
                    clock_time="{:04d} {:02d} {:02d} {:02d} {:02d} {:02d}".format(e[0], e[1], e[2], e[3], e[4], int(e[5]))
                    corr_str = np.array2string(row, separator=', ')
                    with open(file_ssr,'a') as fp:
                        # fp.write(corr_str+"\n")
                        fp.write('> CLOCK {} {}  {:4d} {} \n'.format(clock_time, 0, nsat_clock,"B2B"))
                        for sat_id, data in clock_data.items():
                            fp.write(data)
                    run_clock_update_time=time_clock_sat
                    clock_data={}
                else:
                    # 符合BNC格式的输出
                    str_clock=('{}  {:6d} {:11.4f} {:11.4f} {:11.4f} \n'.format(sat_id, iode,dclk,0.0,0.0))
                    clock_data[sat_id]=str_clock
            if self.subtype == sCSSR.ORBIT:
                iode = self.lc[0].iode[sat]
                dorb = self.lc[0].dorb[sat]  # radial,along-track,cross-track
                if np.isnan(dorb@dorb):
                    continue
                if run_orbit_update_time is None:
                    run_orbit_update_time= self.lc[0].t0[sat][sCType.ORBIT]
                time_orbit_sat=self.lc[0].t0[sat][sCType.ORBIT]
                if timediff(time_orbit_sat,run_orbit_update_time)>1:
                    # epoch update, output the products
                    nsat_orbit=len(orbit_data)
                    e = time2epoch(run_orbit_update_time)
                    orbit_time="{:04d} {:02d} {:02d} {:02d} {:02d} {:02d}".format(e[0], e[1], e[2], e[3], e[4], int(e[5]))
                    corr_str = np.array2string(row, separator=', ')
                    with open(file_ssr,'a') as fp:
                        # fp.write(corr_str+"\n")
                        fp.write('> ORBIT {} {}  {:4d} {} \n'.format(orbit_time, 0, nsat_orbit,"B2B"))
                        for sat_id, data in orbit_data.items():
                            fp.write(data)
                    run_orbit_update_time=time_orbit_sat
                    orbit_data={}
                else:
                    # 符合BNC格式的输出
                    str_orbit=('{}  {:4d} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f}\n'.format(sat_id, iode, dorb[0], dorb[1], dorb[2],0.0,0.0,0.0))
                    orbit_data[sat_id]=str_orbit

    def encode_SP3(self,B2BData0,orb,nav,epoch_time,run_clock_update_time,sp_out,nav_out,file_ssr):
        encodeRTCM=1
        orbit_data = {}
        clock_data = {}
        ns = len(B2BData0.sat_n)
        rs = np.ones((ns, 3))*np.nan
        vs = np.ones((ns, 3))*np.nan
        dts = np.ones((ns, 1))*np.nan

        d_rs = np.ones((ns, 3))*np.nan
        d_dts = np.ones((ns, 1))*np.nan
        peph = peph_t(epoch_time)
        for j, sat in enumerate(B2BData0.sat_n):
            sys, prn = sat2prn(sat)
            sat_id= sat2id(sat)
            if not (sys == uGNSS.GPS or sys==uGNSS.BDS):
                continue
            # cs.iodssr: mask里面解码得到的
            #iodssr_c[sCType.ORBIT] 是轨道部分解码得到的IOD
            if B2BData0.iodssr >= 0 and B2BData0.iodssr_c[sCType.ORBIT] == B2BData0.iodssr:
                # 判断解算的卫星是否在B2b解码的改正数里面
                if sat not in B2BData0.sat_n:
                    continue
            # 判断这颗卫星是否有轨道改正
            if sat not in B2BData0.lc[0].iode.keys():
                continue
            if sat not in B2BData0.lc[0].dorb:
                print("Error")

            iode = B2BData0.lc[0].iode[sat]
            dorb = B2BData0.lc[0].dorb[sat]  # radial,along-track,cross-track
            if B2BData0.lc[0].iodc[sat] == B2BData0.lc[0].iodc_c[sat]:
                dclk = B2BData0.lc[0].dclk[sat]
            else:
                dclk = B2BData0.lc[0].dclk[sat]
                if ~np.isnan(dclk):
                    str_error=("ERROR: Different orbit-clock IOD : "+time2str(epoch_time) +" "+sat_id+" "+\
                               str(B2BData0.lc[0].iodc[sat]) + " " + str(B2BData0.lc[0].iodc_c[sat]))
                    self.log_msg(str_error)
                continue
            if np.isnan(dclk) or np.isnan(dorb@dorb):
                continue

            if sys in B2BData0.nav_mode.keys():
                mode = B2BData0.nav_mode[sys]
            else:
                continue

            eph = findeph(nav.eph, epoch_time, sat, iode=iode, mode=mode)
            if eph is None:
                """
                print("ERROR: cannot find BRDC for {} mode {} iode {} at {}"
                      .format(sat2id(sat), mode, iode, time2str(time)))
                """
                continue

            rs[j, :], vs[j, :], dts[j] = eph2pos(epoch_time, eph, True)
            drel=eph2rel(epoch_time,eph)
            # Along-track, cross-track and radial conversion
            #
            er = vnorm(rs[j, :])
            rc = np.cross(rs[j, :], vs[j, :])
            ec = vnorm(rc)
            ea = np.cross(ec, er)
            A = np.array([er, ea, ec])
            # str_rac="{} {} dorb rac [m] {:14.3f} {:14.3f} {:14.3f} clk [ms] {:12.6f}".format(time2str(run_clock_update_time), sat2id(sat),dorb[0], dorb[1], dorb[2], dclk*1e6)

            # Convert orbit corrections from orbital frame to ECEF
            #
            dorb_e = dorb@A

            # Apply SSR correction
            #
            rs[j, :] -= dorb_e
            dts_brdc=dts[j]*rCST.CLIGHT
            dts[j] += dclk/rCST.CLIGHT  # [m] -> [s]
            # str_xyz="{} {} ssrc xyz [m] {:14.3f} {:14.3f} {:14.3f} clk [ms] {:12.6f}".format(time2str(epoch_time), sat2id(sat),rs[j, 0], rs[j, 1], rs[j, 2], dclk)
            # self.log_msg(str_xyz)

            # 这里的时间其实是最新的钟差
            # dtclk=timediff(epoch_time,B2BData0.lc[0].t0[sat][sCType.CLOCK])
            dtclk=timediff(epoch_time,B2BData0.lc[0].t0[sat][sCType.CLOCK])
            dtorb=timediff(epoch_time,B2BData0.lc[0].t0[sat][sCType.ORBIT])
            str_iod="nav_iod={:4d} mask_iod={:4d} clock_iod={:4d} orbit_iod={:4d}"\
                .format(B2BData0.lc[0].iode[sat],B2BData0.iodssr,B2BData0.lc[0].iodc_c[sat],B2BData0.lc[0].iodc[sat])
            if dtclk>max_clock_delay or dtorb>max_orbit_delay:
                self.log_msg("ERROR: large orbit/clock difference")
                self.log_msg("ERROR Data: "+str_iod)
                B2BData0.deletePRN(sat)
                continue
            else:
                self.log_msg(str_iod)
            peph.pos[sat-1,0]=rs[j,0]
            peph.pos[sat-1,1]=rs[j,1]
            peph.pos[sat-1,2]=rs[j,2]
            peph.pos[sat-1,3]=dts[j,0]-drel
            if sat not in sp_out.sat:
                sp_out.sat.append(sat)

            #encode orbit and clock products
            str_orbit = ('{}  {:4d} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f}\n'. \
                         format(sat_id, iode, dorb[0], dorb[1], dorb[2], 0.0,0.0, 0.0))
            str_clock = ('{}  {:6d} {:11.4f} {:11.4f} {:11.4f} \n'.format(sat_id, iode, dclk, 0.0, 0.0))
            clock_data[sat_id] = str_clock
            orbit_data[sat_id] = str_orbit

        nav_out.peph.append(peph)
        nav_out.ne+=1

        if encodeRTCM==1:
            nsat_clock = len(clock_data)
            nsat_orbit = len(orbit_data)
            e = time2epoch(epoch_time)
            str_time = "{:04d} {:02d} {:02d} {:02d} {:02d} {:02d}".format(e[0], e[1], e[2], e[3], e[4], int(e[5]))
            with open(file_ssr, 'a') as fp:
                fp.write('> CLOCK {} {}  {:4d} {} \n'.format(str_time, 0, nsat_clock, "B2B"))
                for sat_id, data in clock_data.items():
                    fp.write(data)
                fp.write('> ORBIT {} {}  {:4d} {} \n'.format(str_time, 0, nsat_orbit, "B2B"))
                for sat_id, data in orbit_data.items():
                    fp.write(data)
            if nsat_orbit != nsat_clock or nsat_clock<15:
                str_error=">>>>error_time={} nsat_orbit={:4d} nsat_clock={:4d}".format(time2str(epoch_time),nsat_orbit,nsat_clock)
                self.log_msg(str_error)

