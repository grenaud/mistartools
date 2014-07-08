#!/usr/bin/perl


use strict;
use warnings;
use Cwd 'abs_path';

#usage [in mistar file] [population file] [out mistar file]
#where the [population file] is
#
#individual1 pop
#individual2 pop
#...

my %hashPop2Ind;


open(FILE,$ARGV[1]) or die "cannot open ".$ARGV[1];
while(my $line = <FILE>){
  chomp($line);
  #print $line."\n";
  if($line =~ /^(\S+)\s+(\S+)$/){
    print $1."\t".$2."\n";
    if(exists $hashPop2Ind{$2}){
      push(@{$hashPop2Ind{$2}},$1);
    }else{
      $hashPop2Ind{$2} = [$1];
    }
  }else{
    die "wrong line ".$line;
  }
}
close(FILE);

my @path = split("/",abs_path($0));
pop(@path);
my $pathprog = join("/",@path)."/mistarmeld";

my @popFound = (keys %hashPop2Ind);

for(my $keyi=0;$keyi<=$#popFound;$keyi++){
  #foreach my $key (keys %{hashPop2Ind}){
  
  print $pathprog."  ";
  if($keyi == 0){
    print " ".$ARGV[0]." ";
  }else{
    print " /dev/stdin ";
  }
  print " \"".join(",",@{$hashPop2Ind{  $popFound[$keyi] }})."\" \"".$popFound[$keyi]."\" ";
  if($keyi == $#popFound){
    print " |bgzip >  ".$ARGV[2]."\n";
  }else{
    print " | ";
  }
  
}
