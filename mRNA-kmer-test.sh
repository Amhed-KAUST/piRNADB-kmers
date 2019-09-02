##Create working directory
mkdir piRNADBv2
cd piRNADBv2/

###Software to install###
#bedtools ##Not necesary but makes the things easier
#bwa ##Or any other short aligner that produces sam files
#samtools ##Necesary to work with sam/bam files
###

###Starting files####
##Latest C. elegans N2 Worm base annotations##
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS270/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS270.annotations.gff3.gz
##Latest C. elegans genome / Essentially any version since WB235 will be the same##
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS270/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS270.genomic.fa.gz
#Uncompress it
gzip -d c_elegans.PRJNA13758.WS270.genomic.fa.gz
#Index it
bwa index c_elegans.PRJNA13758.WS270.genomic.fa
###

###Command lines###
##Obtain coordinates of all the Wormbase CDS and convert them into BED format. Please note that gff are 1-based format and bed are 0-based##
zcat c_elegans.PRJNA13758.WS270.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t" '{OFS="\t"; if($7 == "-"){str="+"}else{str="-"}; print $1,($4-1),$5,$9,$6,str}' > Celegans-CDS_WB270.bed

###Obtain fasta sequences of CDS in 5` to 3` orientation (strand dependent) ###
bedtools getfasta -s -fi c_elegans.PRJNA13758.WS270.genomic.fa -bed Celegans-CDS_WB270.bed -name+ > Celegans-CDS_WB270.fasta

##Convert each CDS into 12mers and 13mers positioned after the position 2, potentially filtering for uniqueness between the first 14bps and 15 bps disregarding first two positions.
#Please note we believe that 1 to 15 bp alignment is a good measure of strigency
awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-19); i++){print name"separatortoto"substr($0,i,20)"\n"substr($0,i+2,12)}}}' Celegans-CDS_WB270.fasta > Celegans-CDS_WB270.12mers_shift2.fasta
awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-19); i++){print name"separatortoto"substr($0,i,20)"\n"substr($0,i+2,13)}}}' Celegans-CDS_WB270.fasta > Celegans-CDS_WB270.13mers_shift2.fasta

##Look for uniqueness in 12-13 mers with perl script
./QueryInGenomeWithMMinWindow.pl c_elegans.PRJNA13758.WS270.genomic.fa Celegans-CDS_WB270.12mers_shift2.fasta 12 0 1 12 12.tab
./QueryInGenomeWithMMinWindow.pl c_elegans.PRJNA13758.WS270.genomic.fa Celegans-CDS_WB270.13mers_shift2.fasta 13 0 1 13 13.tab

##Convert into fasta files sequences that appear only once, place it as 20mer and remove Ts.
awk -F"\t" '{if($3 == 1){print $1}}' 12.tab | awk -F"separatortoto" '{print $1"\n"$2}' | awk '{if($0 ~ />/){name=$0}else{seq=substr($0,1,3); if(!(seq ~ /T/)){print name"\n"$0}}}' - > Unique-12mers-noT.fasta
awk -F"\t" '{if($3 == 1){print $1}}' 13.tab | awk -F"separatortoto" '{print $1"\n"$2}' | awk '{if($0 ~ />/){name=$0}else{seq=substr($0,1,3); if(!(seq ~ /T/)){print name"\n"$0}}}' - > Unique-13mers-noT.fasta

##Filter for up to two MM in seed region as a 20mer
./QueryInGenomeWithMMinWindow.pl c_elegans.PRJNA13758.WS270.genomic.fa Unique-12mers-noT.fasta 20 2 1 6 Unique-12merShift-noT-2MMseed.tab
./QueryInGenomeWithMMinWindow.pl c_elegans.PRJNA13758.WS270.genomic.fa Unique-13mers-noT.fasta 20 2 1 6 Unique-13merShift-noT-2MMseed.tab

##Convert unique sequences into fasta
awk -F"\t" '{if($3 == 1){print $1"\n"$2}}' Unique-12merShift-noT-2MMseed.tab > Stringent.fasta
awk -F"\t" '{if($3 == 1){print $1"\n"$2}}' Unique-13merShift-noT-2MMseed.tab > LessStringent.fasta

##Map sequences and get bed file
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Stringent.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Stringent.fasta > Stringent.sam
grep "XT:A:U" Stringent.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Stringent.bed

bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa LessStringent.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - LessStringent.fasta > LessStringent.sam
grep "XT:A:U" LessStringent.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > LessStringent.bed



###Please note that any of the fasta files were filtered for GC content; in case that that GC filtering is needed, it can be done in the fasta files with the following command:
#: awk '{if($0 ~ />/){name=$0;}else{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq=$0;k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; GCcontent=(sumC+sumG)/k*100; if((GCcontent >= 30) && (GCcontent <= 45)){print name"_GC-"GCcontent"\n"seq}}}' File.fasta > GCFilt.fasta
## Where 30 is lower boundary and 45 is upper boundary, and File.fasta is input fasta and GCfilt.fasta is name of output fasta.
##Also, this file can be mapped once again with bwa and converted in bed with following command:
#: bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa GCFilt.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - GCFilt.fasta > GCFilt.sam

#: grep "XT:A:U" GCFilt.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > GCFilt.bed

#zcat c_elegans.PRJNA13758.WS270.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t" '{OFS="\t"; if($7 == "-"){str="-"}else{str="+"}; print $1,($4-1),$5,$9,$6,str}' > Celegans-ProperCDS_WB270.bed

#bedtools getfasta -s -fi c_elegans.PRJNA13758.WS270.genomic.fa -bed Celegans-ProperCDS_WB270.bed -name+ > Celegans-ProperCDS_WB270.fasta


zcat c_elegans.PRJNA13758.WS270.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="exon"){print $0}}}' - | grep "Transcript" | awk -F"\t" '{OFS="\t"; if($7 == "-"){str="-"}else{str="+"}; print $1,($4-1),$5,$9":"str":",$6,str}' > Celegans-mRNA_WB270.bed

bedtools getfasta -s -fi c_elegans.PRJNA13758.WS270.genomic.fa -bed Celegans-mRNA_WB270.bed -name > Celegans-mRNA_WB270.fasta

awk '{if($0 ~ />/){name=$0}else{print name"\t"$0}}' Celegans-mRNA_WB270.fasta | awk -F":|\t" '{print $2"\t"$3"\t"$5}' | awk -F"\t" '{if(array[$1] != 0){if($2=="-"){array[$1]=$3""array[$1]}else{array[$1]=array[$1]""$3}}else{array[$1]=$3}}END{for(gen in array){print ">"gen"\n"array[gen]}}' > Celegans-Complete_mRNAs_WB270.fasta


zcat c_elegans.PRJNA13758.WS270.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="mRNA"){print $0}}}' - | awk -F"\t" '{size=($5-$4); print size"\t"$9}' | awk -F"\t" '{print $1";"$2}' | awk -F":|;" '{print $1"\t"$3"\t"$5}' | awk -F"\t" '{if(array[$3]==0){array[$3]=$2; val[$3]=$1}else{if(val[$3]<$1){val[$3]=$1; array[$3]=$2}}}END{for(tra in array){print tra"\t"array[tra]}}' > Gene-longest-Transcript.txt

awk -F"\t" '{print $2}' Gene-longest-Transcript.txt > List-tra.txt

awk '{if($0 ~ />/){name=$0}else{print name"\t"$0}}' Celegans-Complete_mRNAs_WB270.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - List-tra.txt > mRNAReference.fasta


./QueryInGenomeWithMMinWindow2.pl mRNAReference.fasta Stringent.fasta 20 0 1 6 Stringent.tab
./QueryInGenomeWithMMinWindow2.pl mRNAReference.fasta LessStringent.fasta 20 0 1 6 LessStringent.tab

