# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 21:01:49 2021

@author: ruihi
"""

from B2b_HAS_decoder.gnss import *
import numpy as np
from math import pow, sin, cos


NMAX = 10
MAXDTE = 900.0
EXTERR_CLK = 1e-3
EXTERR_EPH = 5e-7


class peph_t:
    def __init__(self, time=None):
        if time is not None:
            self.time = time
        else:
            self.time = gtime_t()
        self.pos = np.ones((uGNSS.MAXSAT, 4))*np.nan
        self.vel = np.ones((uGNSS.MAXSAT, 4))*np.nan
        self.std = np.zeros((uGNSS.MAXSAT, 4))
        self.vst = np.zeros((uGNSS.MAXSAT, 4))


class peph:
    nsat = 0
    nep = 0
    t0 = None
    week0 = -1
    svid = None
    acc = None
    svpos = None
    svclk = None
    svvel = None
    status = 0
    sat=[]
    scl = [0.0, 0.0]

    def __init__(self):
        self.t = None
        self.nmax = 24*12

    def parse_satlist(self, line):
        n = len(line[9:60])//3
        for k in range(n):
            svid = line[9+3*k:12+3*k]
            if int(svid[1:]) > 0:
                self.sat[self.cnt] = id2sat(svid)
                self.cnt += 1

    def parse_acclist(self, line):
        n = len(line[9:60])//3
        for k in range(n):
            acc = int(line[9+3*k:12+3*k])
            if self.cnt < self.nsat:
                self.acc[self.cnt] = acc
            self.cnt += 1

    def parse_sp3(self, fname, nav, opt=0):
        ver_t = ['c', 'd']
        self.status = 0
        v = False
        with open(fname, "r") as fh:
            for line in fh:
                if line[0:3] == 'EOF':  # end of record
                    break
                if line[0:2] == '/*':  # skip comment
                    continue
                if line[0:2] == '* ':  # header of body part
                    self.status = 10

                if self.status == 0:
                    if line[0] != '#':
                        break
                    self.ver = line[1]
                    if self.ver not in ver_t:
                        print("invalid version: {:s}".format(self.ver))
                        break
                    self.flag = line[2]
                    if self.flag not in ('P', 'V'):
                        print("invalid P/V flag: {:s}".format(self.flag))
                        break
                    self.t0 = str2time(line, 3, 27)
                    self.status = 1
                elif self.status == 1:
                    if line[0:2] != '##':
                        break
                    self.week0 = int(line[3:7])
                    # print("week={:4d}".format(self.week0))
                    self.status = 2
                elif self.status == 2:
                    if line[0:2] == '+ ':
                        self.cnt = 0
                        self.nsat = int(line[3:6])
                        self.sat = np.zeros(self.nsat, dtype=int)
                        self.acc = np.zeros(self.nsat, dtype=int)
                        self.parse_satlist(line)
                        self.status = 3
                        continue
                elif self.status == 3:
                    if line[0:2] == '++':
                        self.cnt = 0
                        self.status = 4
                        self.parse_acclist(line)
                        continue
                    self.parse_satlist(line)
                elif self.status == 4:
                    if line[0:2] == '%c':
                        # two %c, two %f, two %i records
                        # Columns 10-12 in the first %c record define
                        # what time system is used for the date/times
                        # in the ephemeris.
                        # c4-5 File type
                        # c10-12 Time System
                        line = fh.readline()  # %c
                        line = fh.readline()  # %f
                        if line[0:2] != '%f':
                            break
                        # Base for Pos/Vel  (mm or 10**-4 mm/sec)
                        # Base for Clk/Rate (psec or 10**-4 psec/sec)
                        self.scl[0] = float(line[3:13])
                        self.scl[1] = float(line[14:26])
                        self.status = 10
                        for _ in range(3):
                            line = fh.readline()
                        continue
                    self.parse_acclist(line)

                if self.status == 10:  # body
                    if line[0:2] == '* ':  # epoch header
                        v = False
                        nav.ne += 1
                        self.cnt = 0
                        peph = peph_t(str2time(line, 3, 27))
                        # ep = time2epoch(peph.time)
                        # print("{:4.0f}/{:02.0f}/{:02.0f} {:2.0f}:{:2.0f}:{:5.2f}"
                        #       .format(ep[0], ep[1], ep[2], ep[3], ep[4], ep[5]))

                        nline = self.nsat if self.flag == 'P' else self.nsat*2
                        for _ in range(nline):
                            line = fh.readline()
                            if line[0] != 'P' and line[0] != 'V':
                                continue

                            svid = line[1:4]
                            sat_ = id2sat(svid)

                            # clk_ev   = line[74] # clock event flag
                            # clk_pred = line[75] # clock pred. flag
                            # mnv_flag = line[78] # maneuver flag
                            # orb_pred = line[79] # orbit pred. flag
                            pred_c = len(line) >= 76 and line[75] == 'P'
                            pred_o = len(line) >= 80 and line[79] == 'P'

                            # x,y,z[km],clock[usec]
                            for j in range(4):
                                if j < 3 and (opt & 1) and pred_o:
                                    continue
                                if j < 3 and (opt & 2) and not pred_o:
                                    continue
                                if j == 3 and (opt & 1) and pred_c:
                                    continue
                                if j == 3 and (opt & 2) and not pred_c:
                                    continue

                                val = float(line[4+j*14:18+j*14])
                                if abs(val-999999.999999) >= 1e-6:
                                    scl = 1e3 if j < 3 else 1e-6
                                    if line[0] == 'P':
                                        v = True
                                        peph.pos[sat_-1, j] = val*scl
                                    elif v:
                                        peph.vel[sat_-1, j] = val*scl*1e-4

                            if len(line) >= 74:
                                for j in range(4):
                                    if j < 3:
                                        slen, scl, ofst = 2, 1e-3, 0
                                    else:
                                        slen, scl, ofst = 3, 1e-12, 1
                                    s = line[61+j*3:61+j*3+slen]
                                    std = int(s) if s[-1] != ' ' else 0
                                    if self.scl[ofst] > 0.0 and std > 0.0:
                                        v = pow(self.scl[ofst], std)*scl
                                        if line[0] == 'P':
                                            peph.std[sat_-1, j] = v
                                        else:
                                            peph.vst[sat_-1, j] = v*1e-4
                    if v:
                        nav.peph.append(peph)

        return nav

    def write_sp3(self, fname, nav):
        """
        Write data to SP3 file
        """

        with open(fname, "w") as fh:

            # Write header section
            #

            # Epoch lines
            #
            t0 = nav.peph[0].time
            e = time2epoch(t0)
            ne = len(nav.peph)

            fh.write("#dP{:04d} {:02d} {:02d} {:02d} {:02d} {:011.8f} {:7d} d+D {:16s}\n"
                     .format(e[0], e[1], e[2], e[3], e[4], e[5], ne, ' '))

            week, secs = time2gpst(t0)
            tstep = timediff(nav.peph[1].time, t0)
            mjd = 44244 + 7*week + int(secs/86400.0)
            fod = time2doy(t0) % 1

            fh.write("## {:04d} {:15.8f} {:14.8f} {:5n} {:15.13f}\n"
                     .format(week, secs, tstep, mjd, fod))

            # Satellite list and accuracy indicators
            #
            self.sat=sorted(self.sat)
            self.nsat=len(self.sat)
            for i in range(int(np.ceil(self.nsat / 17))):

                nsat = "{:4n}".format(self.nsat) if i == 0 else "    "
                prns = [sat2id(s) for s in self.sat[i*17:i*17+17]]
                fh.write('+ {}   {:51s}\n'.format(nsat, ''.join(prns)))

            for i in range(int(self.nsat / 17+1)):

                accs = ['  0' for s in self.sat[i*17:i*17+17]]
                fh.write('++{}   {:51s}\n'.format('    ', ''.join(accs)))

            fh.write(
                '%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n')
            fh.write(
                '%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n')

            fh.write(
                '%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n')
            fh.write(
                '%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n')
            fh.write(
                '%i    0    0    0    0      0      0      0      0         0\n')
            fh.write(
                '%i    0    0    0    0      0      0      0      0         0\n')

            # Comment section
            #
            fh.write('/* \n')

            # Write data section
            #
            for peph in nav.peph:

                e = time2epoch(peph.time)

                fh.write("*  {:04d} {:02d} {:02d} {:02d} {:02d} {:011.8f}\n"
                         .format(e[0], e[1], e[2], e[3], e[4], e[5]))

                for sat in self.sat:

                    if np.isnan(peph.pos[sat-1][0:3]).any():
                        continue

                    clk = 0.999999999999 \
                        if np.isnan(peph.pos[sat-1][3]) else peph.pos[sat-1][3]

                    fh.write("P{:3s} {:13.6f} {:13.6f} {:13.6f} {:13.6f}\n"
                             .format(sat2id(sat),
                                     peph.pos[sat-1][0]*1e-3,
                                     peph.pos[sat-1][1]*1e-3,
                                     peph.pos[sat-1][2]*1e-3,
                                     clk*1e+6))

            # Terminate file
            #
            fh.write('EOF\n')

def Rx(t):
    ct, st = cos(t), sin(t)
    return np.array([[1.0, 0.0, 0.0], [0.0, ct, st], [0.0, -st, ct]])


def Ry(t):
    ct, st = cos(t), sin(t)
    return np.array([[ct, 0.0, -st], [0.0, 1.0, 0.0], [st, 0.0, ct]])


def Rz(t):
    ct, st = cos(t), sin(t)
    return np.array([[ct, st, 0.0], [-st, ct, 0.0], [0.0, 0.0, 1.0]])


def nut_iau1980(t, f):
    nut = np.array([
        [0,   0,   0,   0,   1, -6798.4, -171996, -174.2, 92025,   8.9],
        [0,   0,   2,  -2,   2,   182.6,  -13187,   -1.6,  5736,  -3.1],
        [0,   0,   2,   0,   2,    13.7,   -2274,   -0.2,   977,  -0.5],
        [0,   0,   0,   0,   2, -3399.2,    2062,    0.2,  -895,   0.5],
        [0,  -1,   0,   0,   0,  -365.3,   -1426,    3.4,    54,  -0.1],
        [1,   0,   0,   0,   0,    27.6,     712,    0.1,    -7,   0.0],
        [0,   1,   2,  -2,   2,   121.7,    -517,    1.2,   224,  -0.6],
        [0,   0,   2,   0,   1,    13.6,    -386,   -0.4,   200,   0.0],
        [1,   0,   2,   0,   2,     9.1,    -301,    0.0,   129,  -0.1],
        [0,  -1,   2,  -2,   2,   365.2,     217,   -0.5,   -95,   0.3],
        [-1,   0,   0,   2,   0,    31.8,     158,    0.0,    -1,   0.0],
        [0,   0,   2,  -2,   1,   177.8,     129,    0.1,   -70,   0.0],
        [-1,   0,   2,   0,   2,    27.1,     123,    0.0,   -53,   0.0],
        [1,   0,   0,   0,   1,    27.7,      63,    0.1,   -33,   0.0],
        [0,   0,   0,   2,   0,    14.8,      63,    0.0,    -2,   0.0],
        [-1,   0,   2,   2,   2,     9.6,     -59,    0.0,    26,   0.0],
        [-1,   0,   0,   0,   1,   -27.4,     -58,   -0.1,    32,   0.0],
        [1,   0,   2,   0,   1,     9.1,     -51,    0.0,    27,   0.0],
        [-2,   0,   0,   2,   0,  -205.9,     -48,    0.0,     1,   0.0],
        [-2,   0,   2,   0,   1,  1305.5,      46,    0.0,   -24,   0.0],
        [0,   0,   2,   2,   2,     7.1,     -38,    0.0,    16,   0.0],
        [2,   0,   2,   0,   2,     6.9,     -31,    0.0,    13,   0.0],
        [2,   0,   0,   0,   0,    13.8,      29,    0.0,    -1,   0.0],
        [1,   0,   2,  -2,   2,    23.9,      29,    0.0,   -12,   0.0],
        [0,   0,   2,   0,   0,    13.6,      26,    0.0,    -1,   0.0],
        [0,   0,   2,  -2,   0,   173.3,     -22,    0.0,     0,   0.0],
        [-1,   0,   2,   0,   1,    27.0,      21,    0.0,   -10,   0.0],
        [0,   2,   0,   0,   0,   182.6,      17,   -0.1,     0,   0.0],
        [0,   2,   2,  -2,   2,    91.3,     -16,    0.1,     7,   0.0],
        [-1,   0,   0,   2,   1,    32.0,      16,    0.0,    -8,   0.0],
        [0,   1,   0,   0,   1,   386.0,     -15,    0.0,     9,   0.0],
        [1,   0,   0,  -2,   1,   -31.7,     -13,    0.0,     7,   0.0],
        [0,  -1,   0,   0,   1,  -346.6,     -12,    0.0,     6,   0.0],
        [2,   0,  -2,   0,   0, -1095.2,      11,    0.0,     0,   0.0],
        [-1,   0,   2,   2,   1,     9.5,     -10,    0.0,     5,   0.0],
        [1,   0,   2,   2,   2,     5.6,      -8,    0.0,     3,   0.0],
        [0,  -1,   2,   0,   2,    14.2,      -7,    0.0,     3,   0.0],
        [0,   0,   2,   2,   1,     7.1,      -7,    0.0,     3,   0.0],
        [1,   1,   0,  -2,   0,   -34.8,      -7,    0.0,     0,   0.0],
        [0,   1,   2,   0,   2,    13.2,       7,    0.0,    -3,   0.0],
        [-2,   0,   0,   2,   1,  -199.8,      -6,    0.0,     3,   0.0],
        [0,   0,   0,   2,   1,    14.8,      -6,    0.0,     3,   0.0],
        [2,   0,   2,  -2,   2,    12.8,       6,    0.0,    -3,   0.0],
        [1,   0,   0,   2,   0,     9.6,       6,    0.0,     0,   0.0],
        [1,   0,   2,  -2,   1,    23.9,       6,    0.0,    -3,   0.0],
        [0,   0,   0,  -2,   1,   -14.7,      -5,    0.0,     3,   0.0],
        [0,  -1,   2,  -2,   1,   346.6,      -5,    0.0,     3,   0.0],
        [2,   0,   2,   0,   1,     6.9,      -5,    0.0,     3,   0.0],
        [1,  -1,   0,   0,   0,    29.8,       5,    0.0,     0,   0.0],
        [1,   0,   0,  -1,   0,   411.8,      -4,    0.0,     0,   0.0],
        [0,   0,   0,   1,   0,    29.5,      -4,    0.0,     0,   0.0],
        [0,   1,   0,  -2,   0,   -15.4,      -4,    0.0,     0,   0.0],
        [1,   0,  -2,   0,   0,   -26.9,       4,    0.0,     0,   0.0],
        [2,   0,   0,  -2,   1,   212.3,       4,    0.0,    -2,   0.0],
        [0,   1,   2,  -2,   1,   119.6,       4,    0.0,    -2,   0.0],
        [1,   1,   0,   0,   0,    25.6,      -3,    0.0,     0,   0.0],
        [1,  -1,   0,  -1,   0, -3232.9,      -3,    0.0,     0,   0.0],
        [-1,  -1,   2,   2,   2,     9.8,      -3,    0.0,     1,   0.0],
        [0,  -1,   2,   2,   2,     7.2,      -3,    0.0,     1,   0.0],
        [1,  -1,   2,   0,   2,     9.4,      -3,    0.0,     1,   0.0],
        [3,   0,   2,   0,   2,     5.5,      -3,    0.0,     1,   0.0],
        [-2,   0,   2,   0,   2,  1615.7,      -3,    0.0,     1,   0.0],
        [1,   0,   2,   0,   0,     9.1,       3,    0.0,     0,   0.0],
        [-1,   0,   2,   4,   2,     5.8,      -2,    0.0,     1,   0.0],
        [1,   0,   0,   0,   2,    27.8,      -2,    0.0,     1,   0.0],
        [-1,   0,   2,  -2,   1,   -32.6,      -2,    0.0,     1,   0.0],
        [0,  -2,   2,  -2,   1,  6786.3,      -2,    0.0,     1,   0.0],
        [-2,   0,   0,   0,   1,   -13.7,      -2,    0.0,     1,   0.0],
        [2,   0,   0,   0,   1,    13.8,       2,    0.0,    -1,   0.0],
        [3,   0,   0,   0,   0,     9.2,       2,    0.0,     0,   0.0],
        [1,   1,   2,   0,   2,     8.9,       2,    0.0,    -1,   0.0],
        [0,   0,   2,   1,   2,     9.3,       2,    0.0,    -1,   0.0],
        [1,   0,   0,   2,   1,     9.6,      -1,    0.0,     0,   0.0],
        [1,   0,   2,   2,   1,     5.6,      -1,    0.0,     1,   0.0],
        [1,   1,   0,  -2,   1,   -34.7,      -1,    0.0,     0,   0.0],
        [0,   1,   0,   2,   0,    14.2,      -1,    0.0,     0,   0.0],
        [0,   1,   2,  -2,   0,   117.5,      -1,    0.0,     0,   0.0],
        [0,   1,  -2,   2,   0,  -329.8,      -1,    0.0,     0,   0.0],
        [1,   0,  -2,   2,   0,    23.8,      -1,    0.0,     0,   0.0],
        [1,   0,  -2,  -2,   0,    -9.5,      -1,    0.0,     0,   0.0],
        [1,   0,   2,  -2,   0,    32.8,      -1,    0.0,     0,   0.0],
        [1,   0,   0,  -4,   0,   -10.1,      -1,    0.0,     0,   0.0],
        [2,   0,   0,  -4,   0,   -15.9,      -1,    0.0,     0,   0.0],
        [0,   0,   2,   4,   2,     4.8,      -1,    0.0,     0,   0.0],
        [0,   0,   2,  -1,   2,    25.4,      -1,    0.0,     0,   0.0],
        [-2,   0,   2,   4,   2,     7.3,      -1,    0.0,     1,   0.0],
        [2,   0,   2,   2,   2,     4.7,      -1,    0.0,     0,   0.0],
        [0,  -1,   2,   0,   1,    14.2,      -1,    0.0,     0,   0.0],
        [0,   0,  -2,   0,   1,   -13.6,      -1,    0.0,     0,   0.0],
        [0,   0,   4,  -2,   2,    12.7,       1,    0.0,     0,   0.0],
        [0,   1,   0,   0,   2,   409.2,       1,    0.0,     0,   0.0],
        [1,   1,   2,  -2,   2,    22.5,       1,    0.0,    -1,   0.0],
        [3,   0,   2,  -2,   2,     8.7,       1,    0.0,     0,   0.0],
        [-2,   0,   2,   2,   2,    14.6,       1,    0.0,    -1,   0.0],
        [-1,   0,   0,   0,   2,   -27.3,       1,    0.0,    -1,   0.0],
        [0,   0,  -2,   2,   1,  -169.0,       1,    0.0,     0,   0.0],
        [0,   1,   2,   0,   1,    13.1,       1,    0.0,     0,   0.0],
        [-1,   0,   4,   0,   2,     9.1,       1,    0.0,     0,   0.0],
        [2,   1,   0,  -2,   0,   131.7,       1,    0.0,     0,   0.0],
        [2,   0,   0,   2,   0,     7.1,       1,    0.0,     0,   0.0],
        [2,   0,   2,  -2,   1,    12.8,       1,    0.0,    -1,   0.0],
        [2,   0,  -2,   0,   1,  -943.2,       1,    0.0,     0,   0.0],
        [1,  -1,   0,  -2,   0,   -29.3,       1,    0.0,     0,   0.0],
        [-1,   0,   0,   1,   1,  -388.3,       1,    0.0,     0,   0.0],
        [-1,  -1,   0,   2,   1,    35.0,       1,    0.0,     0,   0.0],
        [0,   1,   0,   1,   0,    27.3,       1,    0.0,     0,   0.0]
    ])

    dpsi = deps = 0.0
    for i in range(106):
        ang = 0.0
        for j in range(5):
            ang += nut[i][j]*f[j]
        dpsi += (nut[i][6]+nut[i][7]*t)*sin(ang)
        deps += (nut[i][8]+nut[i][9]*t)*cos(ang)
    dpsi *= 1E-4*rCST.AS2R  # 0.1 mas -> rad
    deps *= 1E-4*rCST.AS2R
    return dpsi, deps


def time2sec(time):
    ep = time2epoch(time)
    sec = ep[3]*3600.0+ep[4]*60.0+ep[5]
    ep[3] = ep[4] = ep[5] = 0.0
    day = epoch2time(ep)
    return sec, day


def utc2gmst(t: gtime_t, ut1_utc):
    ep2000 = [2000, 1, 1, 12, 0, 0]

    tut = timeadd(t, ut1_utc)
    ut, tut0 = time2sec(tut)
    t1 = timediff(tut0, epoch2time(ep2000))/86400.0/36525.0
    t2 = t1*t1
    t3 = t2*t1
    gmst0 = 24110.54841+8640184.812866*t1+0.093104*t2-6.2E-6*t3
    gmst = gmst0+1.002737909350795*ut

    return np.fmod(gmst, 86400.0)*np.pi/43200.0  # 0 <= gmst <= 2*PI


def orb2ecef(time, rs):
    """
    Rotation matrix from satellite antenna frame to ECEF frame assuming
    standard yaw attitude law
    """

    erpv = np.zeros(5)
    rsun, _, _ = sunmoonpos(gpst2utc(time), erpv, True)
    r = -rs
    ez = r/np.linalg.norm(r)
    r = rsun-rs
    es = r/np.linalg.norm(r)
    r = np.cross(ez, es)
    ey = r/np.linalg.norm(r)
    ex = np.cross(ey, ez)

    return np.array([ex, ey, ez])


def eci2ecef(tutc, erpv):
    ep2000 = [2000, 1, 1, 12, 0, 0]

    tgps = utc2gpst(tutc)
    t = (timediff(tgps, epoch2time(ep2000))+19.0+32.184)/86400.0/36525.0
    t2 = t**2
    t3 = t2*t
    f = ast_args(t)

    ze = (2306.2181*t+0.20188*t2+0.017998*t3)*rCST.AS2R
    th = (2004.3109*t-0.42665*t2-0.041833*t3)*rCST.AS2R
    z = (2306.2181*t+1.09468*t2+0.018203*t3)*rCST.AS2R
    eps = (84381.448-46.8150*t-0.00059*t2+0.001813*t3)*rCST.AS2R
    P = Rz(-z)@Ry(th)@Rz(-ze)

    dpsi, deps = nut_iau1980(t, f)
    N = Rx(-eps-deps)@Rz(-dpsi)@Rx(eps)

    gmst_ = utc2gmst(tutc, erpv[2])
    gast = gmst_+dpsi*np.cos(eps)
    gast += (0.00264*np.sin(f[4])+0.000063*np.sin(2.0*f[4]))*rCST.AS2R

    W = Ry(-erpv[0])@Rx(-erpv[1])
    U = W@Rz(gast)@N@P

    return U, gmst_


def ast_args(t):
    ''' astronomical arguments: f={l,l',F,D,OMG} (rad) '''
    # coefficients for iau 1980 nutation
    fc = np.array([[134.96340251, 1717915923.2178,  31.8792,  0.051635, -0.00024470],
                   [357.52910918,  129596581.0481,  -
                       0.5532,  0.000136, -0.00001149],
                   [93.27209062, 1739527262.8478, -12.7512, -0.001037,  0.00000417],
                   [297.85019547, 1602961601.2090,  -
                       6.3706,  0.006593, -0.00003169],
                   [125.04455501,   -6962890.2665,   7.4722,  0.007702, -0.00005939]])
    f = np.zeros(5)
    tt = np.zeros(4)
    tt[0] = t
    for i in range(1, 4):
        tt[i] = tt[i-1]*t
    for i in range(5):
        f[i] = fc[i][0]*3600.0
        for j in range(4):
            f[i] += fc[i][j+1]*tt[j]
        f[i] = np.fmod(f[i]*rCST.AS2R, 2.0*np.pi)
    return f


def sunmoonpos_eci(tut, rsun=False, rmoon=False):
    ep2000 = [2000, 1, 1, 12, 0, 0]
    t = timediff(tut, epoch2time(ep2000))/86400.0/36525.0
    f = ast_args(t)
    eps = 23.439291-0.0130042*t
    sine, cose = sin(eps*rCST.D2R), cos(eps*rCST.D2R)
    if rsun:
        Ms = 357.5277233+35999.05034*t
        ls = 280.460+36000.770*t+1.914666471 * \
            sin(Ms*rCST.D2R)+0.019994643*sin(2.0*Ms*rCST.D2R)
        rs = rCST.AU*(1.000140612-0.016708617*cos(Ms*rCST.D2R) -
                      0.000139589*cos(2.0*Ms*rCST.D2R))
        sinl, cosl = sin(ls*rCST.D2R), cos(ls*rCST.D2R)
        rsun = rs*np.array([cosl, cose*sinl, sine*sinl])
    if rmoon:
        lm = 218.32+481267.883*t+6.29*sin(f[0])-1.27*sin(f[0]-2.0*f[3]) \
            + 0.66*sin(2.0*f[3])+0.21*sin(2.0*f[0]) - \
            0.19*sin(f[1])-0.11*sin(2.0*f[2])
        pm = 5.13*sin(f[2])+0.28*sin(f[0]+f[2])-0.28 * \
            sin(f[2]-f[0])-0.17*sin(f[2]-2.0*f[3])
        rm = rCST.RE_WGS84/sin((0.9508+0.0518*cos(f[0])+0.0095*cos(f[0]-2.0*f[3])
                               + 0.0078*cos(2.0*f[3])+0.0028*cos(2.0*f[0]))*rCST.D2R)
        sinl, cosl = sin(lm*rCST.D2R), cos(lm*rCST.D2R)
        sinp, cosp = sin(pm*rCST.D2R), cos(pm*rCST.D2R)
        rmoon = rm*np.array([cosp*cosl, cose*cosp*sinl -
                            sine*sinp, sine*cosp*sinl+cose*sinp])
    return rsun, rmoon


def sunmoonpos(tutc, erpv, rsun=False, rmoon=False, gmst=False):
    tut = timeadd(tutc, erpv[2])
    rs, rm = sunmoonpos_eci(tut, rsun, rmoon)
    U, gmst_ = eci2ecef(tutc, erpv)
    if rsun:
        rsun = U@rs
    if rmoon:
        rmoon = U@rm
    if gmst:
        gmst = gmst_
    return rsun, rmoon, gmst