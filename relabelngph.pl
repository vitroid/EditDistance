#!/usr/bin/env perl

# usage: this VMRK < NGPH
# routing2の出力したラベル付け替え情報VMRKを読みこんで、NGPHのラベルを付け替える。

#routing5はサイズの異なる構造でもむりやり照合する。
#その場合
#もし、乱れた構造のほうがノードが少ない場合は、ノードを増やして出力することになる。
#もし、乱れた構造のほうがノードが多い場合は、多いほうのノードで出力する。たぶん、hbroutingでは、多すぎるノードは単に無視される(べき)

USE PYTHON VERSION.

use strict;

open FILE,"<".$ARGV[0];
my @label;
my $nnode; #in ice
while(<FILE>){
  if (/^\@VMRK/){
    $nnode = <FILE>;
    foreach my $i ( 0..$nnode-1 ){
      $_ = <FILE>;
      chomp;
      if ( $_ >= 0 ){
	$label[$_] = $i;
      }
    }
    last;
  }
}

while(<STDIN>){
  if (/^\@NGPH/){
    print;
    my $n = <STDIN>;
    print $nnode;
    while(<STDIN>){
      chomp;
      my ( $x, $y ) = split;
      last if $x < 0;
      $x = $label[$x];
      $y = $label[$y];
      print "$x $y\n";
    }
    print "-1 -1\n";
  }
}
