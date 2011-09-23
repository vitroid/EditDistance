#!/usr/bin/env python

# usage: this VMRK < NGPH
# routing2の出力したラベル付け替え情報VMRKを読みこんで、NGPHのラベルを付け替える。

#routing5はサイズの異なる構造でもむりやり照合する。
#その場合
#もし、乱れた構造のほうがノードが少ない場合は、ノードを増やして出力することになる。
#もし、乱れた構造のほうがノードが多い場合は、多いほうのノードで出力する。たぶん、hbroutingでは、多すぎるノードは単に無視される(べき)

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
