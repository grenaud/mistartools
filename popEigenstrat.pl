#!/usr/bin/perl


use strict;
use warnings;

#usage : popEigenstrat.pl [ind.pop] [pop files]

#this script reads the individual pop file produced
#by mistar2EIGENSTRAT and another file with the following format:
#
#pop=ind1,ind2,...
#
#and produces another individual file where the population
#have been added



my %hashInd2pop;

open(FILE,$ARGV[1]) or die "cannot open ".$ARGV[1];
while(my $line = <FILE>){
  chomp($line);

  if($line =~ /^(.*)=(.*)$/){
    my $pop=$1;
    my $inds=$2;
    #print "inds ".$inds."\n";
    foreach my $ind (split(",",$inds)){
      #print "ind $ind\n";
      if(exists $hashInd2pop{$ind}){
	die "ERROR: The following individual has been found twice: ".$ind." in line $line\n";
      }else{
	$hashInd2pop{$ind} = $pop;
      }
    }

  }else{
    die "ERROR: The following line did not parse ".$line."\n";
  }

}
close(FILE);

open(FILE2,$ARGV[0]) or die "cannot open ".$ARGV[0];
while(my $line = <FILE2>){

  if($line =~ /^(\S+)\s+(\S+)\s+(\S+)$/){
    print $1."\t".$2."\t";
    if(exists $hashInd2pop{$3} ){
      print $hashInd2pop{$3};
    }else{
      print $3;
    }
    print "\n";
  }else{
    die "ERROR: The following line did not parse ".$line."\n";
  }

}
close(FILE2);



