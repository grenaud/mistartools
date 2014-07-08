#!/usr/bin/perl


use strict;
use warnings;
#use Scalar::Util::Numeric qw(isint);
sub thouse{
  my ($num) = @_;
  $num = reverse $num;     # reverse the number's digits
  $num =~ s/(\d{3})/$1,/g; # insert comma every 3 digits, from beginning
  $num = reverse $num;     # Reverse the result
  $num =~ s/^\,//;         # remove leading comma, if any
  return $num;
}



#usage [divergence file] [pop name]



#print 
my $firstLine=1;

my $stringPopname = $ARGV[1];

my @arrayCategories=(  "all",
		       "onlyCpg",
		       "noCpg",
		       "transitions",
		       "transversions",
		       "noDamage"
		    );

my $header = "ind1\tind2\tsites\tAAA\tDDA\tDAA\tADA\tDAA/(DAA+DDA)\tDAA/(DAA+DDA)l\tDAA/(DAA+DDA)h\n";





my $numberOfTargets=0;
open(FILE,$ARGV[0]) or die "cannot open ".$ARGV[0];

while(my $line = <FILE>){
  chomp($line);
  if($line =~ /^(\S+)\-(\S+)\t/){
    #print $2."\n";
    if($2 eq $stringPopname){
      my @array=split("\t",$line);
      my @array2;

      if($firstLine){
	print $header;
	$firstLine=0;
      }


      for(my $i=1;$i<=$#array;$i++){
	if($array[$i] =~ /^\d+$/){
	  #die $array[$i];
	  $array[$i]=thouse($array[$i]);
	}
	push(@array2,$array[$i]);
      }

      my $ind=0;
      for(my $i=0;$i<6;$i++){
	my @array3;

	for(my $j=0;$j<10;$j++){
	  push(@array3,$array2[$ind++]);
	}
	@array3 = splice(@array3,0,7);
	print ($1."\t".$2."\t".$arrayCategories[$i] ."\t". join("\t",@array3)."\n");
      }


    }

     # print $line;
    #}
  }
}
close(FILE);

