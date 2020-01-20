###Pipeline for C. briggase piRNA track based on multiple ones developed for C. elegans ###

###Software ###
#bedtools ##Not necesary but makes the things easier
#bwa ##Or any other short aligner that produces sam files
#samtools ##Necesary to work with sam/bam files
#QueryInGenomeWithMMinWindow.pl #Custom made script
###

##Create working directory
mkdir piRNADB_c_briggase
cd piRNADB_c_briggase


###Starting files####
##Latest C. Briggsae CB4 / Essentially any version since WB254 will be the same##
##Latest C. Briggsae (CB4) Worm base annotations##  ##WS274 at the moment
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS274/species/c_briggsae/PRJNA10731/c_briggsae.PRJNA10731.WS274.annotations.gff3.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS274/species/c_briggsae/PRJNA10731/c_briggsae.PRJNA10731.WS274.genomic.fa.gz

#Uncompress genome
gzip -d c_briggsae.PRJNA10731.WS274.genomic.fa.gz
#Index it
bwa index c_briggsae.PRJNA10731.WS274.genomic.fa
###

###Command lines###
##Obtain coordinates of all the Wormbase CDS and convert them into BED format. Please note that gff are 1-based format and bed are 0-based##
zcat c_briggsae.PRJNA10731.WS274.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="CDS"){print $0}}}' - | awk -F"\t" '{OFS="\t"; if($7 == "-"){str="+"}else{str="-"}; print $1,($4-1),$5,$9,$6,str}' > Cbriggsae-CDS_WB274.bed

###Obtain fasta sequences of CDS in 5` to 3` orientation (strand dependent) ###
bedtools getfasta -s -fi c_briggsae.PRJNA10731.WS274.genomic.fa -bed Cbriggsae-CDS_WB274.bed -name+ > Cbriggsae-CDS_WB274.fasta

##Convert each CDS into n-mers 
#Please note we believe that 2 to 19 bp alignment is a good measure of strigency
awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-19); i++){print name"_"i"separatortoto"substr($0,i,20)"\n"substr($0,i+1,18)}}}' Cbriggsae-CDS_WB274.fasta > CDS-18mers.fasta

##Look for uniqueness in 18 mers with perl script
../QueryInGenomeWithMMinWindow2.pl c_briggsae.PRJNA10731.WS274.genomic.fa CDS-18mers.fasta 18 0 1 18 18.tab

##Convert into fasta files sequences that appear only once, place it as 20mer and remove Ts.
awk -F"\t" '{if($3 == 1){print $1}}' 18.tab | awk -F"separatortoto" '{print $1"\n"$2}' | awk '{if($0 ~ />/){name=$0}else{seq=substr($0,1,3); if(!(seq ~ /T/)){print name"\n"$0}}}' - > Unique-18mers-noT.fasta

##Filter for up to two MM in seed region as a 20mer
../QueryInGenomeWithMMinWindow.pl c_briggsae.PRJNA10731.WS274.genomic.fa Unique-18mers-noT.fasta 20 2 1 6 | awk -F"\t" '{if($3==1){print $1"\n"$2}}' > Unique18mers-2MM-1-6.fasta


###Blast part

##After installing blast, it is required to create a database of the genome to map against
makeblastdb -in c_briggsae.PRJNA10731.WS274.genomic.fa -dbtype 'nucl' -title c_briggsae.PRJNA10731.WS274

##Blast filtered unique sequences using "blast short algorithm"
blastn -db c_briggsae.PRJNA10731.WS274.genomic.fa -query Unique18mers-2MM-1-6.fasta -task "blastn-short" -out Blast18mers-2MM-1-6.tab -outfmt 6 -evalue 1 -num_threads 20

##Get names of seqs that are unique in a set that matches at least 15 bp and star the aligment before the second position
awk -F"\t" '{if((($3*$4)/100)>14){if($7 < 3){array[$1]++}}} END{for (key in array){if(array[key]==1){print key}}}' Blast18mers-2MM-1-6.tab > Blast_names_filt.txt

##Parse output to create a file with the indexes of the unique sequences
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' Unique18mers-2MM-1-6.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - Blast_names_filt.txt > piRNA-seqs.fasta

#Align sequences back to the genome and convert map to a bed file
bwa aln -n 0 -o 0 c_briggsae.PRJNA10731.WS274.genomic.fa piRNA-seqs.fasta | bwa samse c_briggsae.PRJNA10731.WS274.genomic.fa - piRNA-seqs.fasta > piRNA-seqs.sam
grep "XT:A:U" piRNA-seqs.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > piRNA-seqs.bed



###Please note that any of the fasta files were filtered for GC content; in case that that GC filtering is needed, it can be done in the fasta files with the following command:
#: awk '{if($0 ~ />/){name=$0;}else{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq=$0;k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; GCcontent=(sumC+sumG)/k*100; if((GCcontent >= 30) && (GCcontent <= 45)){print name"_GC-"GCcontent"\n"seq}}}' File.fasta > GCFilt.fasta
## Where 30 is lower boundary and 45 is upper boundary, and File.fasta is input fasta and GCfilt.fasta is name of output fasta.
##Also, this file can be mapped once again with bwa and converted in bed with following command:
#: bwa aln -n 0 -o 0 c_briggsae.PRJNA10731.WS274.genomic.fa GCFilt.fasta | bwa samse c_briggsae.PRJNA10731.WS274.genomic.fa - GCFilt.fasta > GCFilt.sam

#: grep "XT:A:U" GCFilt.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > GCFilt.bed




