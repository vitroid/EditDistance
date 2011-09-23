#!/usr/bin/env python
# coding: utf-8

#趣向をがらっと変え、環反転による最適化を行う。
#unsolvedはなくなり、icepenaltyも必要なくなる。
#しかし、収束がさらに遅い。f-ringを反転させるだけだとかなりバリアがでかいに違いない。
#初期配置よりも良い解にまったくいきあたらない(;_;)

# きもちわるいので、引数を逆にする。
# this ice.ngph ice.rngs < distorted.ngph
# ringの数えあげはcountrings3.pyにまかせる。どうせ一回やればいいだけだし。

import sys
import re
import random
import math
import heapq
#backward compatibility with python 2.3
try:
    set()
except NameError:
    from sets import Set as set




def flatten(L):       # Flatten linked list of form [0,[1,[2,[]]]]
    while len(L) > 0:
        yield L[0]
        L = L[1]

#http://code.activestate.com/recipes/119466/
def shortest_path(G, start, end):

    q = [(0, start, ())]  # Heap of (cost, path_head, path_rest).
    visited = set()       # Visited vertices.
    while True:
        (cost, v1, path) = heapq.heappop(q)
        if v1 not in visited:
            visited.add(v1)
            if v1 == end:
                return list(flatten(path))[::-1] + [v1]
            path = (v1, path)
            for (v2, cost2) in G[v1].iteritems():
                if v2 not in visited:
                    heapq.heappush(q, (cost + cost2, v2, path))




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
        #print xyz[0],xyz[1]
        #if not network.has_key(xyz[0]):
        #    network[xyz[0]] = dict()
        #if not network.has_key(xyz[1]):
        #    network[xyz[1]] = dict()
        network[xyz[0]][xyz[1]] = 1
        network[xyz[1]][xyz[0]] = -1




def readAR3A(file):
    line = file.readline()
    #print line,
    n = int(line)
    pattern = re.compile(' +')
    coord = []
    for i in range(0,n):
        line = file.readline()
        xyz = pattern.split(line.strip())
        xyz = map(float,xyz)
        coord.append( xyz );
    return coord


def readBOX3(file):
    line = file.readline()
    #print line,
    pattern = re.compile(' +')
    xyz = pattern.split(line.strip())
    xyz = map(float,xyz)
    return xyz

def saveYaplot( network, ref, box,coord ):
    s = "r 0.1\n"
    for i in network.keys():
        for j in network[i].keys():
            y = 0
            if network[i][j] == 1:
                if ref[i].has_key(j):
                    if ref[i][j] == -1:
                        y = 1
                else:
                    y = 2
                if y > 0:
                    d = list(coord[i])
                    e = [0,0,0]
                    for k in range(0,3):
                        d[k] -= coord[j][k]
                        d[k] -= math.floor( d[k] / box[k] + 0.5 ) * box[k]
                        e[k] = d[k] * 0.2
                        d[k] += coord[j][k]
                        e[k] += coord[j][k]
                    s += "y " + str(y) + "\n"
                    s += "@ " + str(y+3) + "\n"
                    s += " ".join( ["l"] + map(str, coord[j] + d) ) + "\n"
                    s += " ".join( ["c"] + map(str, e) ) + "\n"
    for i in ref.keys():
        for j in ref[i].keys():
            y = 0
            if ref[i][j] == 1:
                if network[i].has_key(j):
                    pass
                else:
                    y = 3
                if y > 0:
                    d = list(coord[i])
                    e = [0,0,0]
                    for k in range(0,3):
                        d[k] -= coord[j][k]
                        d[k] -= math.floor( d[k] / box[k] + 0.5 ) * box[k]
                        e[k] = d[k] * 0.2
                        d[k] += coord[j][k]
                        e[k] += coord[j][k]
                    s += "y " + str(y) + "\n"
                    s += "@ " + str(y+3) + "\n"
                    s += " ".join( ["l"] + map(str, coord[j] + d) ) + "\n"
                    s += " ".join( ["c"] + map(str, e) ) + "\n"
    return s



def readRNGS(file):
    line = file.readline()
    #print line,
    n = int(line)
    pattern = re.compile(' +')
    rings = []
    while True:
        line = file.readline()
        xyz = pattern.split(line.strip())
        #print xyz
        xyz = map(int,xyz)
        if xyz[0] == 0:
            return (n,rings)
        #print xyz[0],xyz[1]
        xyz.remove(xyz[0])
        rings.append( xyz )

def saveNGPH( network ):
    s = "@NGPH\n"
    s += "%d\n" % len(network.keys())
    for i in network.keys():
        for j in network[i].keys():
            if network[i][j] == 1:
                s+= "%s %s\n" % (i,j)
    s += "-1 -1\n"
    return s

#与えられた環が現状、F-ringであるかどうかを判定する。
def fring( ring, network ):
    dir = network[ring[-1]][ring[0]]
    for i in range(1,len(ring)):
        x = ring[i-1]
        y = ring[i]
        if network[x][y] != dir:
            return False
    return True

#辞書にないかもしれない要素をアクセスするための関数
def value(network,x,y):
    if network.has_key(x):
        if network[x].has_key(y):
            return network[x][y]
    return 0


def onestep(beta,rings,network,ref,hamming):
    ring = rings[random.randint(0,len(rings))-1] #-1必要か?問題はないけど
    if fring(ring, network):
        nbond = 0
        penalty = 0
        for i in range(0,len(ring)):
            x=ring[i-1]
            y=ring[i]
            ice = network[x][y]
            dis = value(ref,x,y)
            if dis != 0:
                nbond += 1
                penalty += abs(ice-dis)
        post = nbond*2 - penalty
        #if penalty>0:
        #    print hamming,penalty,post
        if post < penalty or math.exp(-beta*(post-penalty)) > random.random():
            #print beta,hamming,penalty,post
            for i in range(0,len(ring)):
                x=ring[i-1]
                y=ring[i]
                network[x][y] = -network[x][y]
                network[y][x] = -network[y][x]
            hamming += post - penalty
    return hamming
        


#逆行する結合だけをたどって、環が描けるケースがあるはず。
#そのような環は無条件に反転して構わない。
def Rondo( network, ref, hamming ):
	def inversiblepath( path ):
		head = path[-1]
		for i in network[head].keys():
			if network[head][i] == 1:
				if ref[head].has_key(i) and ref[head][i] == -1:
					#go ahead
					for j in range(0,len(path)):
						if path[j] == i:
							#closed ring..
							s = map(str,path)
							sys.stderr.write( " ".join(s) + "\n")
							return path[j:len(path)]
					path.append(i)
					return inversiblepath( path )
		return []
	#sys.stderr.write( "Rondo\n")
	mark = dict()
	for head in range(0,len(network)):
		if not mark.has_key( head ):
			mark[head] = 1
			path = [ head ]
			foundpath = inversiblepath( path )
			if len(foundpath) > 0:
				#print foundpath
				for i in foundpath:
					mark[i] = 1
				#その場で反転してしまう。
				for i in range(0,len(foundpath)):
					x = foundpath[i-1]
					y = foundpath[i]
					network[x][y] = -network[x][y]
					network[y][x] = -network[y][x]
					hamming -= 2
				return hamming
	return hamming



#与えられた有向ネットワークで最も長い経路をさがす。
def pathfinder( network, path ):
    head = path[-1]
    if not network.has_key(head):
        return path
    nexts = network[head].keys()
    if len(nexts) == 0:
        return path
    longestpath = path
    for next in nexts:
        for i in range(0,len(path)):
            if path[i] == next:
                #print "LOOP ALERT",path[i:len(path)] + [next]
                #例外として、探索を中断する。
                return path[i:len(path)] + [next]
        result = pathfinder( network, path + [next] )
        if result[0] == result[-1]:
            #loop exception
            return result
        if len(longestpath) < len(result):
            longestpath = list(result)
    return longestpath


def LineDance( network, ref, hamming, beta ):
    #search ends
    donate = dict()
    accept = dict()
    donate2 = dict()
    for i in network.keys():
        for j in network[i].keys():
            if network[i][j] == 1:
                if not donate2.has_key(i):
                    donate2[i] = dict()
                donate2[i][j] = 1
                if not donate2.has_key(j):
                    donate2[j] = dict()
                donate2[j][i] = 99999
                if ref[i].has_key(j) and ref[i][j] == -1:
                    if not donate.has_key(i):
                        donate[i] = dict()
                    donate[i][j] = 1
                    if not accept.has_key(j):
                        accept[j] = dict()
                    accept[j][i] = 1
    #accept/donateには反転可能な結合のみのネットワークができた。
    heads0 = donate.keys()
    heads = []
    for i in heads0:
        if not accept.has_key(i):
            heads.append(i)
            #print "head:",heads
    if len(heads) > 0:
        path = []
        for head in heads:
            result = pathfinder( donate, [head] )
            if result[0] == result[-1]:
                #loop exception
                result.remove(result[-1])
                #print "LOOP EXCEPTION",result
                for i in range(0,len(result)-1):
                    x = result[i-1]
                    y = result[i]
                    network[x][y] = -network[x][y]
                    network[y][x] = -network[y][x]
                    hamming -= 2
                return hamming
            if len(path) < len(result):
                path = list(result)
        head = path[0]
        tail = path[-1]
        #print path
        #末端をつなぐ最短経路をさがす。Dijkstraを用いる。
        shortcut = shortest_path( donate2, tail,head )
        #print shortcut
        if len(shortcut) <= len(path):
            path += shortcut[1:len(shortcut)-1]
            #print path
            for i in range(0,len(path)):
                x = path[i-1]
                y = path[i]
                lastscore = abs(network[x][y] - value(ref,x,y));
                network[x][y] = -network[x][y]
                network[y][x] = -network[y][x]
                newscore = abs(network[x][y] - value(ref,x,y));
                hamming += newscore - lastscore
            return hamming
    return hamming


#テスト
#とりあえず、全部読みこんでf-ringを仕分けする。
coord = None
box   = None
yaplot=0
if sys.argv[1] == "-y" or sys.argv[1] == "-Y":
    yaplot=1
    if sys.argv[1] == "-Y":
        yaplot=2
    coordfile = open(sys.argv[2], "r")
    sys.argv.remove(sys.argv[0]);
    sys.argv.remove(sys.argv[0]);
    while True:
        line = coordfile.readline()
        if not line:
            break
        if line[0:5] == "@BOX3":
            box = readBOX3(coordfile)
        if line[0:5] == "@AR3A" or line[0:5] == "@NX4A":
            coord = readAR3A(coordfile)
            for i in range(0,len(coord)):
                d = coord[i]
                for k in range(0,3):
                    d[k] -= math.floor( d[k] / box[k] + 0.5 ) * box[k]

icengph = open(sys.argv[1], "r")
icerngs = open(sys.argv[2], "r")
file = sys.stdin
nmol = 0
rings = None
network = None
ref = None
while True:
    line = icengph.readline()
    if not line:
        break
    if line[0:5] == "@NGPH":
        (nmol,network) = readNGPH(icengph)
        break

while True:
    line = icerngs.readline()
    if not line:
        break
    if line[0:5] == "@RNGS":
        (nmol,rings) = readRNGS(icerngs)
        break

while True:
    line = file.readline()
    if not line:
        break
    if line[0:5] == "@NGPH":
        (nmol,ref) = readNGPH(file)
        break

#for ring in rings:
#    print fring( ring, network )

hamming = 0
for x in ref.keys():
    for y in ref[x].keys():
        p = ref[x][y]
        v = value(network,x,y)
        hamming += abs(p-v)
for x in network.keys():
    for y in network[x].keys():
        v = value(ref,x,y)
        if v == 0:
            hamming +=1
hamming /= 2
#print hamming

#random.seed(1)
best = 1000000
loop = 0
lastbest = 0
bestngph = ""
bestyaplot=""
if coord != None:
    bestyaplot = saveYaplot( network, ref, box,coord )
    print bestyaplot
if yaplot==2:
    sys.exit(0)
#    print saveYaplot( network, ref, box,coord )
while True:
    for b in range(30,150):
        beta = b / 100.0
        for j in range(0,1000):
            hamming = onestep(beta,rings,network,ref,hamming)
            if hamming < best:
                sys.stderr.write("%s %s\n" % (beta,hamming))
                best = hamming
                lastbest = loop
                bestngph = saveNGPH( network )
                if coord != None:
                    bestyaplot = saveYaplot( network, ref, box,coord )
                    if yaplot==1:
                        print bestyaplot
#                    print bestyaplot
    hamming = Rondo( network, ref, hamming )
    while True:
        newhamming = LineDance( network, ref, hamming, beta )
        if newhamming == hamming:
            break
        hamming = newhamming
    if hamming < best:
        sys.stderr.write("Rondo/LineDance %s\n" % hamming)
        best = hamming
        lastbest = loop
        bestngph = saveNGPH( network )
        if coord != None and yaplot==1:
            bestyaplot = saveYaplot( network, ref, box,coord )
            if yaplot==1:
                print bestyaplot
#            print bestyaplot
    #sys.stderr.write("#LineDance %s\n" % hamming)
    for b in range(150,30,-1):
        beta = b / 100.0
        for j in range(0,1000):
            hamming = onestep(beta,rings,network,ref,hamming)
            if hamming < best:
                sys.stderr.write("%s %s\n" % (beta,hamming))
                best = hamming
                lastbest = loop
                bestngph = saveNGPH( network )
                if coord != None:
                    bestyaplot = saveYaplot( network, ref, box,coord )
                    if yaplot==1:
                        print bestyaplot
#                if coord != None:
#                    print yaplot( network, ref, box,coord )
    loop += 1
    if best == 0:
        break
    if lastbest + 5 < loop:
        break

if coord == None:
    print "@ETOT\n%s" % best
    print bestngph
elif yaplot==2:
    print bestyaplot
