"""
Galileo HAS correction data decoder

[1] Galileo High Accuracy Service Signal-in-Space
  Interface Control Document (HAS SIS ICD), Issue 1.0, May 2022

"""

import numpy as np
import bitstruct as bs
import galois
from B2b_HAS_decoder.cssrlib import cssr, sCSSR, sCSSRTYPE
from B2b_HAS_decoder.peph import peph,peph_t
from B2b_HAS_decoder.gnss import *
from B2b_HAS_decoder.ephemeris import eph2pos, findeph,eph2rel
from B2b_HAS_decoder.cssrlib import cssr, sCSSR, sCSSRTYPE, sCType

class cssr_has(cssr):
    def __init__(self, foutname=None):  # 初始化属性
        super().__init__(foutname)
        self.MAXNET = 1
        self.cssrmode = sCSSRTYPE.GAL_HAS_SIS
        self.dorb_scl = [0.0025, 0.0080, 0.0080]
        self.dclk_scl = 0.0025
        self.dorb_blen = [13, 12, 12]
        self.dclk_blen = 13
        self.cb_blen = 11
        self.cb_scl = 0.02
        self.pb_blen = 11
        self.pb_scl = 0.01  # cycles

        self.GF = galois.GF(256)

    """
    计算给定n位整数的有符号值
    u: 输入的n位整数值。
    n: 整数的位数。
    scl: 缩放因子，用于对整数值进行缩放
    """

    def sval(self, u, n, scl):
        """ calculate signed value based on n-bit int, lsb """
        invalid = -2 ** (n - 1)
        dnu = 2 ** (n - 1) - 1
        y = np.nan if u == invalid or u == dnu else u * scl
        return y

    '''对于Clock消息提供了所有在Mask块中指定的卫星的时钟更正。这意味着对于Mask中的每个卫星，都会有一个时钟更正值。
     Clock Subset消息则只提供了一个卫星子集的时钟更正。这个子集是通过
     Clock Subset Corrections Block中的一个额外的子集掩码来指定的，
     该掩码指出了哪些卫星在这个子集中。这允许在需要时只更新一部分卫星的时钟更正，而不是全部。'''

    def decode_cssr_clk_sub(self, msg, i=0):
        head, i = self.decode_head(msg, i)  # 调用解码头文件函数解码iodssr
        self.flg_net = False
        if self.iodssr != head['iodssr']:  # 检查消息头部中的iodssr字段是否与当前对象中存储的iodssr字段相匹配如果不匹配，则返回-1表示解码失败。
            return -1

        nclk_sub = bs.unpack_from('u4', msg, i)[0]  # 从输入的msg字节流的位置i开始解包一个4字节的无符号整数，并返回解包钟差子类型。
        i += 4
        for j in range(nclk_sub):
            gnss, dcm = bs.unpack_from('u4u2', msg, i)  # 六个一循环的，前四个赋值给gnss（GPS,Galileo），后二个赋值给dcm表示时钟差倍数
            # （Delta Clock Multipliers），是HAS消息中时钟全集更正块的参数。它指示应用于不同GNSS的时钟差更正的乘数。
            i += 6
            idx_s = np.where(np.array(self.gnss_n) == gnss)  # 从 self.gnss_n 中找到与 gnss 相等的值的索引
            mask_s = bs.unpack_from('u' + str(self.nsat_g[gnss]), msg, i)[0]  # 函数解码从字节流 msg 中提取的卫星掩码
            i += self.nsat_g[gnss]  # 移动到下一个卫星
            idx, nclk = self.decode_mask(mask_s, self.nsat_g[gnss])  # 传递有效的位的索引值，对应卫星数的掩码
            for k in range(nclk):  # 有多少课可用卫星就循环多少次
                dclk = bs.unpack_from('s' + str(self.dclk_blen), msg, i) \
                       * self.dclk_scl * (dcm + 1)  # 解钟差改正数
                self.lc[0].dclk[idx_s[idx[k]]] = dclk  # 将钟差改正数存起来
                i += self.dclk_blen
        return i
      #解码每种信息的头部信息
    def decode_head(self, msg, i, st=-1):

        if st == sCSSR.MASK:
            ui = 0
        else:
            ui = bs.unpack_from('u4', msg, i)[0]  # 解码的第一个元素赋值给ui
            i += 4

        if self.tow0 >= 0:
            self.tow = self.tow0 + self.toh
            if self.week >= 0:
                self.time = gpst2time(self.week, self.tow)  # 计算时间进行时间转换

        head = {'uint': ui, 'mi': 0, 'iodssr': self.iodssr}  # 信息类型 SSR 版本号
        return head, i

    def decode_cssr(self, msg, i=0):
        if self.msgtype != 1:  # only MT=1 defined，判断信息类型，如果不为1解码失败
            print(f"invalid MT={self.msgtype}")
            return False
        # time of hour
        # flags: mask,orbit,clock,clock subset,cbias,pbias,mask_id,iodset_id
        self.toh, flags, res, mask_id, self.iodssr = \
            bs.unpack_from('u12u6u4u5u5', msg, i)  # 解码后的值分别赋值给 时内秒、flags、res、掩码标志 和 iodssr。
        i += 32

        if self.monlevel > 0 and self.fh is not None:
            self.fh.write("##### Galileo HAS SSR: TOH{:6d} flags={:12s} mask_id={:2d} iod_s={:1d}\n"
                          .format(self.toh, bin(flags), mask_id, self.iodssr))

        if self.toh >= 3600:  # 检查toh 是否大于等于3600。如果大于等于3600秒，返回 False 表示解码失败。
            print(f"invalid TOH={self.toh}")
            return False
        '''检查标志位中的第6位，即掩码信息（mask block）是否存在。如果存在，将 mask_id 设置为解码的掩码标识，
         并将 subtype 设置为 sCSSR.MASK。然后调用 decode_cssr_mask 方法来解码掩码部分的HAS信息。'''
        if (flags >> 5) & 1:  #
            self.mask_id = mask_id
            self.subtype = sCSSR.MASK
            i = self.decode_cssr_mask(msg, i)
        '''检查标志位中的第5位，即轨道（orbit block）信息是否存在。如果存在，将 subtype 设置为 sCSSR.ORBIT。
        然后调用 decode_cssr_orb 方法来解码轨道部分的HAS信息获取轨道改正数。如果 monlevel 大于0且文件句柄 fh 不为空，则调用 out_log 方法，输出解码内容。'''
        if (flags >> 4) & 1:  # orbit block
            self.subtype = sCSSR.ORBIT
            i = self.decode_cssr_orb(msg, i)
            if self.monlevel > 0 and self.fh is not None:
                self.out_log()
        '''检查标志位中的第4位，即时钟（clock block）信息是否存在。如果存在，将 mask_id_clk 设置为解码的时钟掩码标识，并将 subtype 设置为 sCSSR.CLOCK。
        然后调用 decode_cssr_clk 方法来解码时钟部分的HAS消息，获取钟差改正数。如果 monlevel 大于0且文件句柄 fh 不为空，则调用 out_log 方法，输出解码内容。'''
        if (flags >> 3) & 1:  # clock block
            self.mask_id_clk = mask_id
            self.subtype = sCSSR.CLOCK
            i = self.decode_cssr_clk(msg, i)
            if self.monlevel > 0 and self.fh is not None:
                self.out_log()
        '''检查标志位中的第3位，即时钟子集块（clock subset block）是否存在。
        如果存在，调用 decode_cssr_clk_sub 方法来解码时钟子集部分的HAS信息。'''
        if (flags >> 2) & 1:  # clock subset block
            i = self.decode_cssr_clk_sub(msg, i)
        '''检查标志位中的第2位，即码偏差（code bias block）是否存在。如果存在，将 subtype 设置为 sCSSR.CBIAS。
        然后调用 decode_cssr_cbias 方法来解码码偏差部分的HAS信息。如果 monlevel 大于0且文件句柄 fh 不为空，则调用 out_log 方法输出解码内容。'''
        if (flags >> 1) & 1:  # code bias block
            self.subtype = sCSSR.CBIAS
            i = self.decode_cssr_cbias(msg, i)
            if self.monlevel > 0 and self.fh is not None:
                self.out_log()
        '''检查标志位中的第1位，即相位偏差（phase bias block）是否存在。如果存在，将 subtype 设置为 sCSSR.PBIAS。
        然后调用 decode_cssr_pbias 方法来解码相位偏差部分的HAS信息。如果 monlevel 大于0且文件句柄 fh 不为空，则调用 out_log 方法输出解码内容。'''
        if (flags >> 0) & 1:  # phase bias block
            self.subtype = sCSSR.PBIAS
            i = self.decode_cssr_pbias(msg, i)
            if self.monlevel > 0 and self.fh is not None:
                self.out_log()

     #提取HAS页面头中的字段
    def decode_has_header(self, buff, i):
        if bs.unpack_from('u24', buff, i)[
            0] == 0xaf3bc3:  # 解码24位的无符号整数，然后与预先定义的值 0xaf3bc3 进行比较。如果相等，表示这是一个空消息，直接返回一组全部为0的值。
            return 0, 0, 0, 0, 0

        hass, res, mt, mid, ms, pid = bs.unpack_from('u2u2u2u5u5u8', buff, i)
        # 如果不是空消息，则继续解包HAS消息的头部信息。这一行使用 unpack_from 方法解包了六个整数，每个整数的位数分别为2、2、2、5、5、8位，
        # 对应于 hass、res、mt、mid、ms 和 pid。
        # 表示HAS状态（HAStatus），是HAS页面头中的2位字段。它指示HAS服务的状态，例如"00"表示测试模式，"01"表示操作模式。
        # mt不为1时解码失败，res是保留信息，：mid表示消息ID（Message ID），是HAS页面头中的5位字段。它唯一标识正在传输的消息，ms（Message Size）确实表示未编码消息的大小，以页面数表示
        # pid表示页面ID（Page ID），是HAS页面头中的8位字段。它唯一标识正在传输的HAS编码页面。

        ms += 1
        return hass, mt, mid, ms, pid

    # 这个方法非常重要，因为它可以从接收的这些页面中恢复出完整的HAS消息，进而提取出卫星更正信息。
    def decode_has_page(self, idx, has_pages, gMat, ms):
        """ HPVRS decoding for RS(255,32,224) """
        HASmsg = bytes()
        k = len(idx)
        '''使用生成矩阵和页面内容对HAS进行解码。首先，通过索引列表从 has_pages 中取出相应的页面内容，
        然后将这个页面内容表示为Galois域上的矩阵 Wd，尺寸为k x 53。同时，从生成矩阵中取出与索引列表对应的部分，然后对该部分进行逆矩阵运算，
        得到逆矩阵 Dinv，尺寸为 k x k。接着，将 Dinv 乘以 Wd 得到解码后的消息 Md，尺寸为 k x 53。最后，将 Md 转换为字节对象并赋值给 HASmsg。'''
        if k >= ms:
            Wd = self.GF(has_pages[idx, :])  # kx53
            Dinv = np.linalg.inv(self.GF(gMat[idx, :k]))  # kxk
            Md = Dinv @ Wd  # decoded message (kx53)
            HASmsg = np.array(Md).tobytes()

        return HASmsg
    def log_msg(self,msg):
        if self.monlevel > 0 and self.fh is not None:
            self.fh.write(msg+"\n")

    def encode_SP3(self,HASData0,orb,nav,epoch_time,sp_out,nav_out,file_ssr):
        max_orbit_delay = 300
        max_clock_delay = 30
        encodeRTCM=1
        orbit_data = {}
        clock_data = {}
        ns = len(HASData0.sat_n)
        rs = np.ones((ns, 3))*np.nan
        vs = np.ones((ns, 3))*np.nan
        dts = np.ones((ns, 1))*np.nan

        d_rs = np.ones((ns, 3))*np.nan
        d_dts = np.ones((ns, 1))*np.nan
        peph = peph_t(epoch_time)
        for j, sat in enumerate(HASData0.sat_n):
            sys, prn = sat2prn(sat)
            sat_id= sat2id(sat)
            if not (sys == uGNSS.GPS or sys==uGNSS.GAL):
                continue

            if HASData0.iodssr_c[sCType.ORBIT] == HASData0.iodssr:
                if sat not in HASData0.sat_n:
                    continue
            else:
                continue

            iode = HASData0.lc[0].iode[sat]
            dorb = HASData0.lc[0].dorb[sat]  # radial,along-track,cross-track

            if HASData0.cssrmode == sCSSRTYPE.GAL_HAS_SIS:  # HAS only
                if HASData0.mask_id != HASData0.mask_id_clk:  # mask has changed
                    if sat not in HASData0.sat_n_p: #clock sat list
                        continue
            else:
                print("=======>Error: unrecongnized format for PPP")

            dclk = HASData0.lc[0].dclk[sat]

            if np.isnan(dclk) or np.isnan(dorb@dorb):
                continue

            if sys in HASData0.nav_mode.keys():
                mode = HASData0.nav_mode[sys]
            else:
                continue

            eph = findeph(nav.eph, epoch_time, sat, iode=iode, mode=mode)
            if eph is None:
                """
                print("ERROR: cannot find BRDC for {} mode {} iode {} at {}"
                      .format(sat2id(sat), mode, iode, time2str(time)))
                """
                continue
            if HASData0.mask_id != HASData0.mask_id_clk:  # mask has changed
                self.log_msg("ERROR: not matching id for orbit and clock, oribt_id= "+str(HASData0.mask_id)+" clock_id= "+str(HASData0.mask_id_clk)+"  "+time2str(epoch_time))
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
            # Convert orbit corrections from orbital frame to ECEF
            #
            dorb_e = dorb@A

            # Apply SSR correction
            #
            rs[j, :] -= dorb_e
            dts_brdc=dts[j]*rCST.CLIGHT
            dts[j] += dclk/rCST.CLIGHT  # [m] -> [s]


            peph.pos[sat-1,0]=rs[j,0]
            peph.pos[sat-1,1]=rs[j,1]
            peph.pos[sat-1,2]=rs[j,2]
            peph.pos[sat-1,3]=dts[j,0]-drel
            if sat not in sp_out.sat:
                sp_out.sat.append(sat)


            dtclk=timediff(epoch_time,HASData0.lc[0].t0[sat][sCType.CLOCK])
            dtorb=timediff(epoch_time,HASData0.lc[0].t0[sat][sCType.ORBIT])
            str_iod="nav_iod={:4d} mask_iod={:4d} ".format(HASData0.lc[0].iode[sat],HASData0.iodssr)
            if dtclk>max_clock_delay or dtorb>max_orbit_delay:
                self.log_msg("ERROR: large orbit/clock difference")
                self.log_msg("ERROR Data: "+str_iod)
                HASData0.deletePRN(sat)
                continue

            #encode orbit and clock products
            str_orbit = ('{}  {:4d} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f}\n'. \
                         format(sat_id, iode, dorb[0], dorb[1], dorb[2], 0.0,0.0, 0.0))
            str_clock = ('{}  {:6d} {:11.4f} {:11.4f} {:11.4f} \n'.format(sat_id, iode, dclk, 0.0, 0.0))
            clock_data[sat_id] = str_clock
            orbit_data[sat_id] = str_orbit

        nsat_clock = len(clock_data)
        nsat_orbit = len(orbit_data)
        if nsat_orbit != nsat_clock or nsat_clock<15:
            str_error=">>>>error_time={} nsat_orbit={:4d} nsat_clock={:4d}".format(time2str(epoch_time),nsat_orbit,nsat_clock)
            self.log_msg(str_error)
            return

        nav_out.peph.append(peph)
        nav_out.ne+=1

        if encodeRTCM==1:
            e = time2epoch(epoch_time)
            str_time = "{:04d} {:02d} {:02d} {:02d} {:02d} {:02d}".format(e[0], e[1], e[2], e[3], e[4], int(e[5]))
            with open(file_ssr, 'a') as fp:
                fp.write('> CLOCK {} {}  {:4d} {} \n'.format(str_time, 0, nsat_clock, "HAS"))
                for sat_id, data in clock_data.items():
                    fp.write(data)
                fp.write('> ORBIT {} {}  {:4d} {} \n'.format(str_time, 0, nsat_orbit, "HAS"))
                for sat_id, data in orbit_data.items():
                    fp.write(data)
