#!/usr/bin/env perl
#hamming$B5wN%$K$h$kJ?6QE*$J3H;6$NMM;R$r8+$k!"$"$k$$$O5wN%9TNs$r7W;;$9$k!#(B
#$BG[Ns$N;H$$J}$rJQ99$9$k!#CY$9$.$k$N$G!#(B
#$BJ?@.(B13$BG/(B7$B7n(B11$BF|(B($B?e(B)$B%P%0H/8+(B
use strict;

use vars qw(%first %second);

sub distance{
    my($broken)=@_;
    my($tot1,$tot2,$common,$i);
    if($broken){
	#first$B$NJ}$r4p=`$K$9$k!#(B
	@_=keys %second;
	@_=keys %first;
	$tot1=$#_+1;
	foreach $i ( @_ ){
	    if(defined $second{$i}){
		$common++;
	    }
	}
	if($tot1==0){
	    return 0;
	}
	return $common/$tot1;
    }else{
	@_=keys %second;
	$tot2=$#_+1;
	foreach $i ( @_ ){
	    if(defined $first{$i}){
	    }
	    else{
	      print STDERR $i,"\n";
	    }
	}
	print STDERR "\n";
	@_=keys %first;
	$tot1=$#_+1;
	foreach $i ( @_ ){
	    if(defined $second{$i}){
		$common++;
	    }
	    else{
	      print STDERR $i,"\n";
	    }
	}
	return $tot1+$tot2-2*$common;
    }
}

#-u$B$,M?$($i$l$?$iL58~%0%i%U$H$_$J$9(B
my $u=0;
#-d$B$,M?$($i$l$?$iA4G[CV4V$N5wN%9TNs$r=PNO$9$k!#(B
my $dm=0;
#-p$B$,M?$($i$l$?$i5wN%9TNs$r(Bpgm$B7A<0$G=PNO$9$k!#(B
my $pgm=0;
#-r n$B$,M?$($i$l$?$i!"(Bn$BHVL\$N%9%J%C%W%7%g%C%H$r86E@$H$7$?5wN%$r5a$a$k!#(B
my $ref=-1;
#$BG\N($,M?$($i$l$?$i(Bpgm$B7A<0$NG;EY$rD4@a$9$k!#(B
my $mul=1;
#-i$B$,M?$($i$l$?$i!"3F%9%F%C%W4V$N5wN%$N$_$r=PNO$9$k!#(B
my $int=0;
#-b$B$,M?$($i$l$?$i!"@ZCG$5$l$?(BHB$B$N$_$G5wN%$rDj5A$9$k!#(B
my $broken=0;

while($_=shift){
    if($_ eq "-u"){
	$u=1;
    }elsif($_ eq "-b"){
	$broken=1;
    }elsif($_ eq "-d"){
	$dm=1;
    }elsif($_ eq "-i"){
	$int=1;
	$dm=0;
	$pgm=0;
	$ref=-1;
    }elsif($_ eq "-p"){
	$pgm=1;
	$dm=1;
    }elsif($_ eq "-r"){
	$ref=shift;
    }else{
	$mul=$_;
    }
}

my @g=();
while(<STDIN>){
    if(/\@NGPH/){
	my $n=<STDIN>;
	my %s;
	undef %s;
	while(1){
	    $_=<STDIN>;
	    chomp;
	    last if /-/;
	    s/ +/ /g;
	    s/^\s//;
	    $s{$_}=1;
	    #print STDERR "[",$_, "]\n";
	    if($u){
		chomp;
		$s{join(" ",reverse split)}=1;
	    }
	}
	push(@g,\%s);
    }
}

my @dist;
my @count;

if($dm&&$pgm){
    printf "P2\n%d %d\n%d\n",$#g+1,$#g+1,255;
}
undef %first;
if($ref>=0){
    %first=%{$g[$ref]};
    for(my $i=0;$i<=$#g;$i++){
	%second=%{$g[$i]};
	my $d=distance($broken);
	if($u){
	    $d/=2;
	}
	print "$ref $i $d\n";
    }
}elsif($int){
    for(my $i=0;$i<=$#g;$i++){
	%second=%{$g[$i]};
	if($i>0){
	    $ref=$i-1;
	    my $d=distance($broken);
	    if($u){
		$d/=2;
	    }
	    print "$ref $i $d\n";
	}
	%first=%second;
    }
}else{
    for(my $i=0;$i<=$#g;$i++){
	undef %first;
	my $ii;
	%first=%{$g[$i]};
	#print STDERR "[$i]";
	if($pgm){
	    for(my $j=0;$j<$i;$j++){
		print "0 ";
	    }
	}
	for(my $j=$i;$j<=$#g;$j++){
	    %second=%{$g[$j]};
	    my $d=distance($broken);
	    if($u){
		$d/=2;
	    }
	    my $dd=int($d*$mul);
	    if($pgm){
		if($dd>255){
		    $dd=255;
		}
		print "$dd ";
	    }else{
		print "$i $j $d\n";
	    }
	    my $k=$j-$i;
	    $count[$k]++;
	    $dist[$k]+=$d;
	}
	print "\n";
    }
}

if(!($dm||$int||$ref)){
    for(my $i=0;$i<$#count;$i++){
	print $i," ",$dist[$i]/$count[$i],"\n";
    }
}
