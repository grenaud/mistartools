#!/usr/bin/perl


use strict;
use warnings;
use Data::Dumper;

warn "reading chromosomes\n";
my $cmd = "zcat ".$ARGV[0]." |grep -v \"^#\" |cut -f 1 |uniq";
my @allchr = split("\n",`$cmd`);

foreach my $chr  (@allchr){
  warn "Processing chr ".$chr."\n";
  my $cmd2 = "tabix -h ".$ARGV[0]." ".$chr." |bgzip -c > ".$ARGV[0].".".$chr.".gz\n";
  #print $cmd2;
  my $out = `$cmd2`;
}




