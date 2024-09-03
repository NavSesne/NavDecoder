#!/usr/bin/env python3
#
#  uinit test for sdr_ldpc.py
#
import sys, time
sys.path.append('../python')
import numpy as np
from B2b_HAS_decoder import sdr_ldpc
from B2b_HAS_decoder import sdr_nb_ldpc
import bitstruct as bs
from binascii import unhexlify,hexlify

# Convert the LDPC(162,91)
def decode_LDPC(vi_row):
    if isinstance(vi_row['nav'], np.ndarray):
        buff=vi_row['nav'][0]
    else:
        buff=vi_row['nav']
    buff=buff.decode('utf-8')
    buff =  buff[:-2]
    # buff = [item.decode('utf-8') for item in buff]
    # buff =  buff[0][:-2]
    length=len(buff)
    # buff = buff.decode('utf-8')
    SF = read_hex(buff)
    SF=SF[12:]
    ok, ng = 0, 0
    t = time.time()

    err_data = SF.copy()

    dec_data, nerr = sdr_ldpc.decode_LDPC('BCNV3', err_data)

    if np.all(SF[:486] == dec_data):
        #print('test_08 (%3d) OK: nerr=%3d' % (n, nerr))
        ok += 1
    else:
        #print('test_08 (%3d) NG: nerr=%3d' % (n, nerr))
        ng += 1
    hex_txt=hex_str(dec_data)
    if len(hex_txt) % 2 == 1:
        hex_txt += '0'
    # hex_txt1=hexlify(dec_data).decode('utf-8')
    buff = unhexlify(hex_txt)
    return buff
# pack bits to uint8 ndarray ---------------------------------------------------
def pack_bits(data, nz=0):
    if nz > 0:
        data = np.hstack([[0] * nz, data])
    N = len(data)
    buff = np.zeros((N + 7) // 8, dtype='uint8')
    for i in range(N):
        buff[i // 8] |= (data[i] << (7 - i % 8))
    return buff

# read HEX strings -------------------------------------------------------------
def read_hex(str0):
    N = len(str0) * 4
    data = np.zeros(N, dtype='uint8')
    for i in range(N):
        data[i] = (int(str0[i // 4], 16) >> (3 - i % 4)) & 1
    return data

# data to HEX strings ----------------------------------------------------------
def hex_str(data):
    str = ''
    hex = 0
    for i in range(len(data)):
        hex = (hex << 1) + data[i]
        if i % 4 == 3:
            str += '%1X' % (hex)
            hex = 0
    return str
