#!/usr/bin/perl


use strict;
use warnings;

#usage gphocs file [fraction to keep]

my $fraction = $ARGV[1];
open(FILE,$ARGV[0]) or die "cannot open ".$ARGV[0];
my $line = <FILE>;
#print $line;#do not print till later
my $stringToPrint="";

if($line =~ /^\d+$/){
  # number of records
}else{
  $stringToPrint.="\n";
}

my $keptLocus=0;
my $state=0; #0 before block, 1=in block
my $lenblock;
my $lastghost=0;
my $flushedlastghost=0;
my $keepLocus=0;

while($line = <FILE>){
  chomp($line);

  if($state==0){
    if(length($line) == 0){
      if($keepLocus){
	#print $line."\n";
	$stringToPrint.=$line."\n";
      }

    }else{
      if($line =~ /locus(\d+)\s(\d+)\s(\d+)$/){

	if(rand()<=$fraction){
	  $keepLocus=1;
	  $keptLocus+=1;
	}else{
	  $keepLocus=0;
	}

	if($keepLocus){
	  #print "locus".$1." ".($2+1)." ".$3."\n";
	  $stringToPrint.= "locus".$keptLocus." ".($2)." ".$3."\n";
	}

	$state    =1;
	$lastghost=0;
	$lenblock =$3;
	$flushedlastghost=0;
      }
    }
  } else {
    if ($state==1) {
      if (length($line) == 0) {	#end of state
	$flushedlastghost=1;
	$state=0;
	if ($keepLocus) {
	  #print $line."\n";
	  $stringToPrint.= $line."\n";
	}
      }else{
	if ($keepLocus) {
	  #print $line."\n";
	  $stringToPrint.= $line."\n";
	}
      }
    }else{
      die "internal error\n";
    }
  }

}
close(FILE);

if(!$flushedlastghost){
  #print "ghost".($lastghost+1)."\t"."N"x$lenblock."\n";
  $stringToPrint.= "ghost".($lastghost+1)."\t"."N"x$lenblock."\n";
}

print $keptLocus."\n\n".$stringToPrint."\n";
