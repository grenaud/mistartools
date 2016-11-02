#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

my $mock=0;

sub runcmd{
  my ($cmdtorun) = @_;

  warn  "running cmd ". $cmdtorun."\n";
  if($mock != 1){
    my @argstorun = ( "bash", "-c", $cmdtorun );

    if(system(@argstorun) != 0){
      die "system  cmd $cmdtorun failed: $?"
    }else{
    }
  }

}



sub fileExistsExec{
  my ($exeFile) = @_;

  if (!( -x $exeFile)) {
    die "Executable ".$exeFile." does not exist\n";
  }
}

sub fileExists{
  my ($exeFile) = @_;

  if (!( -e $exeFile)) {
    die "Executable ".$exeFile." does not exist\n";
  }
}


my $lengthmerge =0;
my $minlength   =1000;
my $outputf   ="none";

my $help;

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\n
 This script is a wrapper to convert mistar files into the seq format
 of g-phocs. Normally, you have to run mistar2gphocs using a bed file
 of regions to consider. This script finds these regions.

\n\n usage:\t".$0." <options> [mst file] \n\n".

" Options:\n".

 "\t--mock\t\t\t\t\tDo nothing, just print the commands that will be run\n".
  "\t-l [X bp]\t\t\t\tConsider overlap between regions if they are Xbp away (default ".$lengthmerge." bp)\n".
  "\t-m [X bp]\t\t\t\tOnly print regions longer than Xbp away (default ".$minlength." bp)\n".
  "\t-o [output]\t\t\t\tOutput file (default: ".$outputf.")\n".

    "\n\n";
  die;
}




 usage() if ( @ARGV < 1 or
             ! GetOptions('help|?' => \$help, 'mock' => \$mock, 'l=i' => \$lengthmerge, 'm=i' => \$minlength,'o=s' => \$outputf)
          or defined $help );

if($outputf eq "none"){
  die "ERROR: output file must be specified";
}

my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);

my $mistar2bed        = $pathdir."/mistar2bed";
my $mistar2gphocs     = $pathdir."/mistar2gphocs";
my $inputMSTfile      = $ARGV[ $#ARGV  ];

fileExistsExec($mistar2bed);
fileExistsExec($mistar2gphocs);
fileExists($inputMSTfile);

my $cmd;

#call mistar2bed





$cmd=$mistar2bed." ".$inputMSTfile." ";

#detect mergeBed
my $cmdmergebeddetect="which mergeBed";
my $mergebebpath=`$cmdmergebeddetect`;
chomp($mergebebpath);
if(!$mergebebpath){
  die "ERROR: cannot detect command mergeBed, please install it and make it accessible in your path\n";
}


if($lengthmerge>0){
  $cmd = $cmd ." | ".$mergebebpath." -d ".$lengthmerge." -i /dev/stdin "
}
#print $lengthmerge;

#filter by length
if($minlength>0){
  $cmd = $cmd ." | awk '{if( (\$3-\$2)>=".$minlength." ){print \$0}}' ";
}

my $outputftemp=$outputf."_";

#call mistar2gphocs
$cmd = $cmd ." | ". $mistar2gphocs ."  ".$inputMSTfile." /dev/stdin > ".$outputftemp;
runcmd($cmd);


my $numberofLoci=0;
open(FILE,$outputftemp) or die "cannot open ".$outputftemp."\n";
while(my $line = <FILE>){
  chomp($line);
  if($line =~ /^locus(\d+)/){
    $numberofLoci++;
  }
}
close(FILE);



open(FILEO,">".$outputf) or die "cannot write to ".$outputf."\n";
print FILEO $numberofLoci."\n\n";

open(FILE,$outputftemp) or die "cannot open ".$outputftemp."\n";

while(my $line = <FILE>){
  print FILEO $line;
}
close(FILE);

close(FILEO);

unlink($outputftemp) or die "cannot remove temp file ".$outputftemp."\n";


warn "mistar2gphocsWrapper.pl finished successfully, wrote data to ".$outputf."\n";
