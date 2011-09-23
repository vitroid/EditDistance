#!/usr/bin/env python

# usage: this VMRK < NGPH
# routing2�ν��Ϥ�����٥��դ��ؤ�����VMRK���ɤߤ���ǡ�NGPH�Υ�٥���դ��ؤ��롣

#routing5�ϥ������ΰۤʤ빽¤�Ǥ�����ȹ礹�롣
#���ξ��
#�⤷����줿��¤�Τۤ����Ρ��ɤ����ʤ����ϡ��Ρ��ɤ����䤷�ƽ��Ϥ��뤳�Ȥˤʤ롣
#�⤷����줿��¤�Τۤ����Ρ��ɤ�¿�����ϡ�¿���ۤ��ΥΡ��ɤǽ��Ϥ��롣���֤�hbrouting�Ǥϡ�¿������Ρ��ɤ�ñ��̵�뤵���(�٤�)

import sys
import re

def readNGPH(file):
    line = file.readline()
    #print line,
    n = int(line)
    pattern = re.compile(' +')
    network = dict()
    for i in range(n):
        network[i] = dict()
    while True:
        line = file.readline()
        xyz = pattern.split(line.strip())
        #print xyz
        xyz = map(int,xyz)
        if xyz[0] < 0:
            return (n,network)
        network[xyz[0]][xyz[1]] = 1
        network[xyz[1]][xyz[0]] = -1


def readVMRK(file):
    line = file.readline()
    #print line,
    n = int(line)
    vmrk = []
    for i in range(0,n):
        line = file.readline()
        vmrk.append(int(line));
    return vmrk


def invert(list):
    tsil = [0] * len(list)
    for i in range(len(list)):
        if list[i] >= 0:
            tsil[list[i]] = i
    return tsil


def seekandread( file, tag, reader ):
    while True:
        line = file.readline()
        if line == "":
            return None,None
        result = tag.match( line )
        if result:
            return reader(file)

            



fileVMRK = open(sys.argv[1], "r")
tagVMRK = re.compile("^@VMRK")
tagNGPH = re.compile("^@NGPH")

while True:
    (n, network) = seekandread( sys.stdin, tagNGPH, readNGPH )
    if not n:
        break
    vmrk = seekandread( fileVMRK,  tagVMRK, readVMRK )
    if not vmrk:
        break
    label = invert( vmrk )
    print "@NGPH"
    print len(vmrk)
    for x,line in network.iteritems():
        for y,value in line.iteritems():
            if value > 0:
                print label[x],label[y]
    print -1,-1
