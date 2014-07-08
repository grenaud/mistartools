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
print "\\documentclass[a4paper,10pt]{article}\n";
print "\\usepackage[landscape]{geometry}\n";
print "\\usepackage[dvips]{graphicx,epsfig}\n";
print "\\usepackage{color}\n";
print " \n";
print "\n";
print "\\begin{document}\n";



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

my $header =
"\\begin{tabular}{".
"|l|l|l|lllllll".
"}\n".
"ind1 & ind2 & sites & AAA & DDA & DAA & ADA & DAA        & DAA        & DAA         \\\\ \n".
"     &      &       &     &     &     &     & (DAA+DDA)  & (DAA+DDA)l & (DAA+DDA)h  \\\\ \n \\hline \n";
#"ind1 & ind2 & sites & AAA & DDA & DAA & ADA & DAA        & DAA        & DAA        & ADA       & ADA        & ADA        \\\\ \n".
#"     &      &       & AAA & DDA & DAA & ADA & (DAA+DDA)  & (DAA+DDA)l & (DAA+DDA)h & (ADA+DDA) & (ADA+DDA)l & (ADA+DDA)h \\\\ \n \\hline \n";


;

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
	#print "\\begin{tabular}{";
	#print "|l|l|l|llllllllll";
	#print "}\n";
	#print "ind1 & ind2 & sites & AAA & DDA & DAA & ADA & DAA/(DAA+DDA)& DAA/(DAA+DDA)l & DAA/(DAA+DDA)h & ADA/(ADA+DDA)& ADA/(ADA+DDA)l & ADA/(ADA+DDA)h \\\\";
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
	print ($1." & ".$2." & ".$arrayCategories[$i] ." & ". join(" & ",@array3)." \\\\ \n");
      }


      #my $string= join(" & ",@array2);
      #chomp($string);
      #print ($1." & ".$2." & ".$string." \\\\ \n");
      if($numberOfTargets == 4){
	print "\\end{tabular}\n";
	print "\\newpage\n";

	print $header;
	#print "\\begin{tabular}{";
	#print "|l|l|l|"."l"x($#array);
	#print "|l|l|l|llllllllll";
	#print "}\n";
	$numberOfTargets =0;
      }else{
	$numberOfTargets++;
      }
      print "\\hline \n";
    }

     # print $line;
    #}
  }
}
close(FILE);


print "\\end{tabular}\n";
print "\\end{document}\n";
