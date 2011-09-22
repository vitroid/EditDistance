#icc$B$K$7$F$_$F$b(B1$B3dDxEY$7$+B.$/$J$i$J$+$C$?(B
#CXX=icpc
CXX=g++
CXXFLAGS=-static -O3 -static #-tpp7 -xW -parallel -gcc-name=g++-3.4 -gcc-version=340
#CXXFLAGS=-O3 -static -tpp7 -xW -parallel #-gcc-name=gcc-3.4 -gcc-version=340
LDFLAGS=-pthread
#DMTX$B7A<0$GA4MWAG$r=PNO!#(B4096$BJ,;R$@$H%U%!%$%k$,5pBg$K$J$k!#(B
#CXXDEBUGFLAGS=-g
all: routing2++.x       
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS) 
%.x: %.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) 
test: ice.ar3a defect.ar3a
	./routing3 -y ice.ar3a < defect.ar3a > converge.yap
test2: int.all
%.all:
	make $*.route2++.log $*.route2++.vmrk $*.route2++.ngph $*.route2++.yap $*.hbroute8+.ngph $*.hbroute8+.yap
#routing2$B$K$h$j!"=E?40LCV$,6aIU$/$h$&$K%i%Y%k$r8r49$9$k!#(B
%.route2.log: ice.ar3a %.ar3a
	./routing2+ -l -5 ice.ar3a < $*.ar3a > $@
%.route2++.log: ice.ar3a %.ar3a
	./routing2++ -l -5 ice.ar3a < $*.ar3a > $@
#$B=E?40\F07PO)$N2D;k2=(B
%.route2++.yap: ice.ar3a %.ar3a %.route2++.ar3a
	cat $*.route2++.ar3a | ./routing2++ -l -5 -y ice.ar3a > $@
%.route2++.Yap: ice.ar3a %.ar3a %.route2++.ar3a
	cat $*.route2++.ar3a | ./routing2++ -l -5 -Y ice.ar3a > $@
#$B:GE,$J%i%Y%j%s%0$r<h$j=P$9!#(B
%.route2.vmrk: %.route2.log
	sed -n -e '1,/#####/d' -e '/VMRK/,$$p' $< > $@
%.route2++.vmrk: %.route2++.log
	sed -n -e '1,/#####/d' -e '/VMRK/,$$p' $< > $@
%.ar3a: %.log
	sed -e '1,/#####/d' -e '/VMRK/,$$d' $< > $@
#$B$=$l$K9g$o$;$F7k9g$r(Brelabel$B$9$k!#(B
%.route2.ngph: %.route2.vmrk %.hb
	perl ./relabelngph.pl $*.route2.vmrk < $*.hb > $@
%.route2.ngph: %.route2.vmrk %.ngph
	perl ./relabelngph.pl $*.route2.vmrk < $*.ngph > $@
%.route2++.ngph: %.route2++.vmrk %.ngph
	perl ./relabelngph.pl $*.route2++.vmrk < $*.ngph > $@
#relabel$B8e$N%0%i%U$O!"L58~%0%i%U$H$7$F$OI9$N%0%i%U$H$+$J$j6a$$$O$:!#(B
#$B$=$3$G!"L58~%0%i%U$H$7$F>H9g$7$F5wN%$r5a$a$k!#(B
%.uhd: ice.hb %.ngph
	cat $^ | ./hamming.pl -u | sed -n -e 2,2p 
#$BM-8~%0%i%U$G$N5wN%$r5a$a$k!#$3$3$G$O!"(B%.route2.ngph$B$N7k9g>pJs$rB:=E$7!"$=$l$K$"$o$;$F(Bice.hb$B$N$[$&$N7k9g$r:GE,2=$9$k!#(B2$B<T$N%M%C%H%o!<%/$N%O%_%s%05wN%$r$G$-$k$@$1>.$5$/$7$D$D!"(Bice rule$B$r$_$?$9$h$&$K$9$k!#(B
%.hbroute6.ngph: %.route2.ngph ice.hb
	perl ./hbrouting6.pl $*.route2.ngph < ice.hb > $@
%.hbroute8.ngph: %.route2.ngph ice.hb
	python ./hbrouting8.py ice.hb ice.rngs < $*.route2.ngph > $@
#%.hbroute8+.ngph: %.route2.ngph ice.hb
#	python ./hbrouting8+.py ice.hb ice.rngs < $*.route2.ngph > $@
%.hbroute8+.ngph: %.route2++.ngph ice.hb
	python ./hbrouting8+.py ice.hb ice.rngs < $*.route2++.ngph > $@
#%.hbroute8+.yap: %.route2.ngph ice.hb
#	python ./hbrouting8+.py -y ice.ar3a ice.hb ice.rngs < $*.route2.ngph > $@
%.hbroute8+.yap: %.hbroute8+.ngph %.route2++.ngph
	python ./hbrouting8+.py -y ice.ar3a $*.hbroute8+.ngph ice.rngs < $*.route2++.ngph > $@
%.hbroute8+.Yap: %.hbroute8+.ngph %.route2++.ngph
	python ./hbrouting8+.py -Y ice.ar3a $*.hbroute8+.ngph ice.rngs < $*.route2++.ngph > $@
%.hbroute8+.yel: %.hbroute8+.ngph %.route2++.ngph
	python ./yellowchain.py $*.hbroute8+.ngph < $*.route2++.ngph > $@
%.route2++.bla: %.ar3a %.route2++.vmrk blackchain
	cat $*.route2++.vmrk $*.ar3a | ./blackchain ice.ar3a > $@
#$B7k9gAH$_49$(7PO)$N2D;k2=(B
#layer 3($BNP(B) $B7k9gDI2C(B
#layer 4($B2+(B) $BJQ2=$J$7(B
#layer 5($B@V(B) $B7k9g@ZCG(B
#layer 6($B@D(B) $B7k9gH?E>(B
%.hbroute6.yap: %.route2.ngph %.route2.ar3a %.hbroute6.ngph t+n22.pl
	cat $*.route2.ngph $*.route2.ar3a | perl ./t+n22.pl -bc BondDiff -d $*.hbroute6.ngph @NGPH -pp BondDiff -d $*.hbroute6.ngph @NGPH -bl BondDiff -d $*.hbroute6.ngph @NGPH  -o OutputYaplot -p PaletteYaplot0 > $@

%: %.in
	pwd=`pwd`; sed -e "s|__INCPATH__|$$pwd|" $< > $@
%.rngs: %.ngph
	../CountRings/countrings2 -C 8 < $< > $@
%.edit: %.route2++.ngph %.hbroute8+.ngph
	python ./editdistance.py $^ > $@
#mark displaced molecules
%.corr.vmrk: ice.ar3a %.route2++.ar3a
	./correspond ice.ar3a < $*.route2++.ar3a > $@
%.corr.yap: %.route2++.ngph %.route2++.ar3a %.corr.vmrk
	~/src/Tools/TrajTools/Offset 0 -6 14 <  $*.route2++.ar3a | cat $*.route2++.ngph -  |  perl ~/src/Tools/Traj2Mol3/t+n22.pl -bl BondLayerVMRK $*.corr.vmrk @VMRK > $@

