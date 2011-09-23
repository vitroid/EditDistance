#!/usr/bin/env perl

# usage: this VMRK < NGPH
# routing2�ν��Ϥ�����٥��դ��ؤ�����VMRK���ɤߤ���ǡ�NGPH�Υ�٥���դ��ؤ��롣

#routing5�ϥ������ΰۤʤ빽¤�Ǥ�����ȹ礹�롣
#���ξ��
#�⤷����줿��¤�Τۤ����Ρ��ɤ����ʤ����ϡ��Ρ��ɤ����䤷�ƽ��Ϥ��뤳�Ȥˤʤ롣
#�⤷����줿��¤�Τۤ����Ρ��ɤ�¿�����ϡ�¿���ۤ��ΥΡ��ɤǽ��Ϥ��롣���֤�hbrouting�Ǥϡ�¿������Ρ��ɤ�ñ��̵�뤵���(�٤�)

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
