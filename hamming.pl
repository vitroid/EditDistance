#!/usr/bin/env perl
#hamming距離による平均的な拡散の様子を見る、あるいは距離行列を計算する。
#配列の使い方を変更する。遅すぎるので。
#平成13年7月11日(水)バグ発見
use strict;

use vars qw(%first %second);

sub distance{
    my($broken)=@_;
    my($tot1,$tot2,$common,$i);
    if($broken){
	#firstの方を基準にする。
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

#-uが与えられたら無向グラフとみなす
my $u=0;
#-dが与えられたら全配置間の距離行列を出力する。
my $dm=0;
#-pが与えられたら距離行列をpgm形式で出力する。
my $pgm=0;
#-r nが与えられたら、n番目のスナップショットを原点とした距離を求める。
my $ref=-1;
#倍率が与えられたらpgm形式の濃度を調節する。
my $mul=1;
#-iが与えられたら、各ステップ間の距離のみを出力する。
my $int=0;
#-bが与えられたら、切断されたHBのみで距離を定義する。
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
