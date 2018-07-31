#!/usr/bin/env python

import sys
from os import path

if (not path.isfile(sys.argv[1])):
    print('Usage: t-j_input.py SYSTEM_FILE.')
    sys.exit()

inp_com = open(sys.argv[1], 'r')
inp_hei = open('model.heis.def', 'w')
inp_hub = open('model.hubb.def', 'w')

def do_pipe(dst_fid):
    while True:
        inp_fln = inp_com.readline()
        if inp_fln.strip() == '&END':
            break
        dst_fid.write(inp_fln.split('#')[0])
    dst_fid.write('\n\n')

def do_pipe2(dst_fid1, dst_fid2):
    while True:
        inp_fln = inp_com.readline()
        if inp_fln.strip() == '&END':
            break
        dst_fid1.write(inp_fln.split('#')[0])
        dst_fid2.write(inp_fln.split('#')[0])
    dst_fid1.write('\n\n')
    dst_fid2.write('\n\n')

def do_pipe0():
    while True:
        inp_fln = inp_com.readline()
        if inp_fln.strip() == '&END':
            break

while True:
    card_title = inp_com.readline().strip()
    if card_title == '':
        break
    if card_title.strip() == '&COMM' or card_title.strip() == '&COMMON':
        do_pipe2(inp_hei, inp_hub)
    elif card_title.strip() == '&HUB' or card_title.strip() == '&HUBBARD':
        do_pipe(inp_hub)
    elif card_title.strip() == '&HEIS' or card_title.strip() == '&HEISENBERG':
        do_pipe(inp_hei)
    elif card_title.strip()[0] == '&':
        do_pipe0()
        print('Warning: Card ' + card_title + ' is not ignored.')

inp_com.close()
inp_hei.close()
inp_hub.close()

