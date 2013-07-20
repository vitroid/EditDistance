# -*- coding: utf-8 -*-
#CXX=icpc
CXX=g++
CXXFLAGS=-O3
LDFLAGS=    #-pthread
all: routing2++.x
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS) 
%.x: %.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) 
test: int.all
%.all:
	make -k $*.route2++.log $*.route2++.vmrk $*.route2++.ngph $*.route2++.yap $*.hbroute8+.ngph $*.hbroute8+.yap
#Estimate the geometrical edit distance (GED) by routing2++.x
%.route2++.log: samples/ice.ar3a samples/%.ar3a routing2++.x
	./routing2++.x -l -5 samples/ice.ar3a < samples/$*.ar3a > $@
#Extract the optimized path info from log file of routing2++.x
%.route2++.vmrk: %.route2++.log
	sed -n -e '/VMRK/,$$p' $< > $@
#Relabel the hydrogen bonds in order that the total displacements be the smallest according to the result by routing2++.x
%.route2++.ngph: %.route2++.vmrk samples/%.hb
	perl ./relabelngph.py $*.route2++.vmrk < samples/$*.hb > $@
#Visualize the rgeometrical edit path for YaPlot
%.route2++.yap: samples/ice.ar3a %.route2++.ar3a
	cat $*.route2++.ar3a | ./routing2++.x -l -5 -y samples/ice.ar3a > $@
%.route2++.Yap: samples/ice.ar3a %.ar3a %.route2++.ar3a
	cat $*.route2++.ar3a | ./routing2++ -l -5 -Y samples/ice.ar3a > $@
#Estimate the topological edit distance (TED) by hbrouting8+.py
%.hbroute8+.ngph: %.route2++.ngph samples/ice.hb
	python ./hbrouting8+.py samples/ice.hb samples/ice.rngs < $*.route2++.ngph > $@
#Visualize the topological edit path for YaPlot
%.hbroute8+.yap: %.hbroute8+.ngph %.route2++.ngph
	python ./hbrouting8+.py -y samples/ice.ar3a $*.hbroute8+.ngph samples/ice.rngs < $*.route2+




#Obsolete / Unused #####################################################
%.ar3a: %.log
	sed -n -e '/#####/,/VMRK/p' $< | sed -e 's/^@VMRK//' > $@
#それに合わせて結合をrelabelする。
%.route2.ngph: %.route2.vmrk %.hb
	perl ./relabelngph.pl $*.route2.vmrk < $*.hb > $@
%.route2.ngph: %.route2.vmrk %.ngph
	perl ./relabelngph.pl $*.route2.vmrk < $*.ngph > $@
#relabel後のグラフは、無向グラフとしては氷のグラフとかなり近いはず。
#そこで、無向グラフとして照合して距離を求める。
%.uhd: samples/ice.hb %.ngph
	cat $^ | ./hamming.pl -u | sed -n -e 2,2p 
#有向グラフでの距離を求める。ここでは、%.route2.ngphの結合情報を尊重し、それにあわせてsamples/ice.hbのほうの結合を最適化する。2者のネットワークのハミング距離をできるだけ小さくしつつ、ice ruleをみたすようにする。
%.hbroute6.ngph: %.route2.ngph samples/ice.hb
	perl ./hbrouting6.pl $*.route2.ngph < samples/ice.hb > $@
%.hbroute8.ngph: %.route2.ngph samples/ice.hb
	python ./hbrouting8.py samples/ice.hb samples/ice.rngs < $*.route2.ngph > $@
#%.hbroute8+.ngph: %.route2.ngph samples/ice.hb
#	python ./hbrouting8+.py samples/ice.hb samples/ice.rngs < $*.route2.ngph > $@
%.hbroute8+.ngph: %.route2++.ngph samples/ice.hb
	python ./hbrouting8+.py samples/ice.hb samples/ice.rngs < $*.route2++.ngph > $@
#%.hbroute8+.yap: %.route2.ngph samples/ice.hb
#	python ./hbrouting8+.py -y samples/ice.ar3a samples/ice.hb samples/ice.rngs < $*.route2.ngph > $@
%.hbroute8+.yap: %.hbroute8+.ngph %.route2++.ngph
	python ./hbrouting8+.py -y samples/ice.ar3a $*.hbroute8+.ngph samples/ice.rngs < $*.route2++.ngph > $@
%.hbroute8+.Yap: %.hbroute8+.ngph %.route2++.ngph
	python ./hbrouting8+.py -Y samples/ice.ar3a $*.hbroute8+.ngph samples/ice.rngs < $*.route2++.ngph > $@
%.hbroute8+.yel: %.hbroute8+.ngph %.route2++.ngph
	python ./yellowchain.py $*.hbroute8+.ngph < $*.route2++.ngph > $@
%.route2++.bla: %.ar3a %.route2++.vmrk blackchain
	cat $*.route2++.vmrk $*.ar3a | ./blackchain samples/ice.ar3a > $@
#結合組み換え経路の可視化
#layer 3(緑) 結合追加
#layer 4(黄) 変化なし
#layer 5(赤) 結合切断
#layer 6(青) 結合反転
%.hbroute6.yap: %.route2.ngph %.route2.ar3a %.hbroute6.ngph t+n22.pl
	cat $*.route2.ngph $*.route2.ar3a | perl ./t+n22.pl -bc BondDiff -d $*.hbroute6.ngph @NGPH -pp BondDiff -d $*.hbroute6.ngph @NGPH -bl BondDiff -d $*.hbroute6.ngph @NGPH  -o OutputYaplot -p PaletteYaplot0 > $@

%: %.in
	pwd=`pwd`; sed -e "s|__INCPATH__|$$pwd|" $< > $@
%.rngs: %.ngph
	../CountRings/countrings2 -C 8 < $< > $@
%.edit: %.route2++.ngph %.hbroute8+.ngph
	python ./editdistance.py $^ > $@
#mark displaced molecules
%.corr.vmrk: samples/ice.ar3a %.route2++.ar3a
	./correspond samples/ice.ar3a < $*.route2++.ar3a > $@
%.corr.yap: %.route2++.ngph %.route2++.ar3a %.corr.vmrk
	~/src/Tools/TrajTools/Offset 0 -6 14 <  $*.route2++.ar3a | cat $*.route2++.ngph -  |  perl ~/src/Tools/Traj2Mol3/t+n22.pl -bl BondLayerVMRK $*.corr.vmrk @VMRK > $@

