#!/usr/bin/env python
# coding: utf-8

#2つのネットワークを読みこんで、切断、結合、反転の個数を出力するのみ。

import sys
import re


#辞書にないかもしれない要素をアクセスするための関数
def value(network,x,y):
    if network[x].has_key(y):
        return network[x][y]
    return 0




def readNGPH(file):
    line = file.readline()
    #print line,
    n = int(line)
    pattern = re.compile(' +')
    network = dict()
    while True:
        line = file.readline()
        xyz = pattern.split(line.strip())
        #print xyz
        xyz = map(int,xyz)
        if xyz[0] < 0:
            return (n,network)
        #print xyz[0],xyz[1]
        if not network.has_key(xyz[0]):
            network[xyz[0]] = dict()
        if not network.has_key(xyz[1]):
            network[xyz[1]] = dict()
        network[xyz[0]][xyz[1]] = 1
        network[xyz[1]][xyz[0]] = -1




icengph = open(sys.argv[1], "r")
refngph = open(sys.argv[2], "r")
ice = None
ref = None
while True:
    line = icengph.readline()
    if not line:
        break
    if line[0:5] == "@NGPH":
        (nmol,ice) = readNGPH(icengph)
        break

while True:
    line = refngph.readline()
    if not line:
        break
    if line[0:5] == "@NGPH":
        (nmol,ref) = readNGPH(refngph)
        break

cut=0
connect=0
flip=0

for x in ref.keys():
    for y in ref[x].keys():
        p = ref[x][y]
        v = value(ice,x,y)
        if p == -v:
            flip+= 1
        elif p == v:
            pass
        elif v == 0:
            connect += 1
for x in ice.keys():
    for y in ice[x].keys():
        v = value(ref,x,y)
        if v == 0:
            cut += 1

#flipはflipの必要数の2倍=edit距離に直してある。
print cut/2, connect/2, flip
