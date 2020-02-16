#!/usr/bin/perl

#####Documentation######
$num_args = $#ARGV + 1;
if ($num_args != 6) {
  print "\nUsage: QueryInGenomeWithMMinWindow.pl Ref_fasta Query_fasta k-mer_length mm_in_window Pos-start_window Pos-end_window\n\n";
  print "\tRef_fasta: Fasta file of genomic sequence to map. It can be multi-fasta\n";
  print "\tQuery_seq: Set of k-mers that are going to be mapped\n";
  print "\tk-mer_length: Size of subsequences that the reference file will be partitioned into. The k-mer length has to be the same size as the length of each query seq\n";
  print "\tmm_in_window: Number of missmatches allowed in within a defined region of the query\n";
  print "\tPos-start_window: Position that the window starts. The first bp is considered to be in the position 1\n";
  print "\tPos-end_window: Position that the window ends.\n";

  print "\n\nNote: This program relies on storing data in RAM and can be too greedy for long k-mers (e.g. 4**20 subsequences can be stored for a k-mer of length 20)\n";
  exit;
}
############

##Reading and verifying input commands
my $Ref=$ARGV[0];
my $Query=$ARGV[1];
my $k=$ARGV[2];
my $mm=$ARGV[3];
my $pwin1=$ARGV[4];
my $pwin2=$ARGV[5];

###If the following files are not present -> stop the program
if (!(-e $Ref."")){die $Ref." file not found!\n";}
if (!(-e $Query."")){die $Name." file not found!\n";}
if ($pwin2 > $k ){die "Window size cannot be larger than the k-mer length\n";}
if ($mm > $k ){die "Total number or missmatches cannot be larger than the k-mer length\n";}
if ($pwin1 < 1 ){die "The window start position has to be positive\n";}
if ($pwin2 < 1 ){die "The window end position has to be positive\n";}
if ($pwin2 < $pwin1 ){die "The window end position has to be larger than the start position\n";}

$pwin1 = $pwin1 - 1;
$pwin2 = $pwin2 - 1; 
#@dict=("A","T","C","G");

##Reading and hashing genomic sequence
my $string="";
open(op,$Ref) or die "cannot open Ref file\n";
while($line=<op>){
chomp($line);
if(($line =~ m/^>/)){$string=$string."__________"; next;}
$string=$string.$line
}
close(op);

$Genomic = $string;
$Revcomp = reverse $Genomic;
$Revcomp =~ tr/ACGTacgt/TGCAtgca/;

#Hash both strands
#print "\nIndexing genome\n";
my %hash;
for($i=0;$i<=(length($Genomic)-$k);$i++){
$hash{substr($Genomic,$i,$k)}++;
$hash{substr($Revcomp,$i,$k)}++;
}

##Read each query and perform actions
my $Qseq;
my @array;
my $sarray;
my $val;
open(op,$Query) or die "cannot open query file\n";
#print "\nReading sequences\n";
##my $con=0;
while($line=<op>){
##$con++;
chomp($line);
if($line =~ m/^>/){$Qseq=$line}else{
if(length($line) != $k){die $Qseq." query sequence has a different size (".length($Qseq).") than k-mer-length (".$k.").\n";}
#if($hash{$line}>1){next;}
#$val=$hash{$line};
#my $flag=1;
#print $Qseq."\t".$line."\t".$val."\n"}}

##If mm larger than 0 go for array approach otherwise print hashed value
$val = 0;
if($mm == 0){
  $val=$hash{$line};
  }else{
    #print "Line to test".$line."\n";
    @array = produceSeqs($line,$pwin1,$pwin2,$mm);
    foreach $finito(@array){
      #print "\nHola:".$finito."\n";
      $val = $val + $hash{$finito}; 
      if(length($finito) != $k){die "Amhed you did something weird with the sequence ".$line.' producing '.$finito}}
  }

print $Qseq."\t".$line."\t".$val."\n";
#if(flag ==1){print $Qseq."\t".$line."\t".$val."\n";}
} #Stop if not header
##if(!($con % 5)){print "\n".$con." sequences analyzed\n";}
} ##Stop reading file
close(op);

exit;

###Sub function: Produce $line $start $end $number of mismatches
sub produceSeqs {
  #print "\nA\n";
  @seq=split('',$_[0]);
  #print @seq;
  #print "\n";
  $left="";
  $right="";
  if($_[1] > 0){$left= $left.join('',@seq[0..($_[1]-1)])}
  if(($_[2]+1) < scalar(@seq)){$right= $right.join('',@seq[($_[2]+1)..(scalar(@seq)-1)])}
  #print "Left:".$left."\n";
  #print "Right:".$right."\n";
  @temp = repseq($_[3],@seq[$_[1]..$_[2]]);
  @seqs=();
  foreach $tem (@temp){push(@seqs,$left.$tem.$right)}
  #print "\nSeqsA\n";
  #print @seqs;
  return(@seqs);
}

sub repseq{
  #print "\nB Args\n";
  @args=@_;
  #print @args;
  $rep=$args[0];
  #print "\nRepetitions:".$rep;
  @letters=@args[1..(scalar(@args)-1)];
  #print "\n Args again\n";
  #print @args;
  #print "\nSize:".scalar(@args);
  #print "\nMaxIndex:".$#args;
  #print "\nletters:\n";
  #print @letters;
  #@dict=("A","T","C","G");
  @res=();
  push(@res,join('',@letters));
  my %rest;
  #$size = keys %rest;
  #print "\nAvantResB\n";
  #print $size;
  while($rep != 0){
    @res=perseq(@res);
    $rep--;
  }
  foreach $s (@res){$rest{$s}++}
  #print "\nResB\n";
  #$size = keys %rest;
  #print $size;
  return(keys(%rest));
}

sub perseq{
  #print "\nC\n";
  @fin=();
  foreach $sequence (@_){
    #print "\nSequence inside C:\n";
    #print "-".$sequence."-";
    @pat=split('',$sequence);
    for($j=0; $j< scalar(@pat);$j++){
      @tmep=@pat;
      $tmep[$j]="A";
      push(@fin,join("",@tmep));
      $tmep[$j]="T";
      push(@fin,join("",@tmep));
      $tmep[$j]="C";
      push(@fin,join("",@tmep));
      $tmep[$j]="G";
      push(@fin,join("",@tmep));
    }
  }
  return(@fin);
}
