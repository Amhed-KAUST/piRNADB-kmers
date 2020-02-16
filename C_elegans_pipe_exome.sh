###Pipeline for C. briggase piRNA track based on multiple ones developed for C. elegans ###

###Software ###
#bedtools ##Not necesary but makes the things easier
#bwa ##Or any other short aligner that produces sam files
#samtools ##Necesary to work with sam/bam files
#QueryInGenomeWithMMinWindow.pl #Custom made script
###

mkdir piRNADB_c_elegans_exome
cd piRNADB_c_elegans_exome

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

zcat c_elegans.PRJNA13758.WS270.annotations.gff3.gz | awk -F"\t" '{if($2=="WormBase"){if($3=="mRNA"){print $0}}}' - | awk -F"\t" '{OFS="\t"; if($7 == "-"){str="-"}else{str="+"}; print $1,($4-1),$5,$9,$6,str}' > Celegans-mRNA_WB270.bed


###Obtain fasta sequences of CDS in 5` to 3` orientation (strand dependent) ###
bedtools getfasta -s -fi c_elegans.PRJNA13758.WS270.genomic.fa -bed Celegans-CDS_WB270.bed -name+ > Celegans-CDS_WB270.fasta

bedtools merge -i Celegans-mRNA_WB270.bed -s -c 6 -o distinct > tmp

mv tmp Celegans-mRNA_WB270.bed

bedtools getfasta -s -fi c_elegans.PRJNA13758.WS270.genomic.fa -bed Celegans-mRNA_WB270.bed -name+ > Celegans-mRNA_WB270.fasta

### Unique sites only
##Convert each CDS into 20-mers 
#0 MM
awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-19); i++){print name"_"i"separatortoto"substr($0,i,20)"\n"substr($0,i,20)}}}' Celegans-CDS_WB270.fasta > CDS-20mers.fasta

##Look for uniqueness in 20 mers with perl script (hash method)
../QueryInGenomeWithMMinWindow2.pl Celegans-mRNA_WB270.fasta CDS-20mers.fasta 20 0 1 20 20.tab

##Convert into fasta files sequences that appear only once, place it as 20mer and remove Ts.
awk -F"\t" '{if($3 == 1){print $1}}' 20.tab | awk -F"separatortoto" '{print $1"\n"$2}' | awk '{if($0 ~ />/){name=$0}else{seq=substr($0,1,3); if(!(seq ~ /T/)){print name"\n"$0}}}' - > Unique-20mers-noT.fasta

##Obtain positions of uniquely mapped by bwa
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Unique-20mers-noT.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Unique-20mers-noT.fasta > Cel_exo_0MM_piRNA-seqs.sam
grep "XT:A:U" Cel_exo_0MM_piRNA-seqs.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_exo_0MM_piRNA-seqs.bed

###1 mismatch by perl strategy only
##Convert into 18-mers then add 1 MM
#1 MM
awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-19); i++){print name"_"i"separatortoto"substr($0,i,20)"\n"substr($0,i+1,18)}}}' Unique-20mers-noT.fasta > CDS-18mers.fasta

##Look for uniqueness in 18 mers with perl script (hash method)
../QueryInGenomeWithMMinWindow2.pl Celegans-mRNA_WB270.fasta CDS-18mers.fasta 18 0 1 18 18.tab

##Convert into fasta files sequences that appear only once, place it as 20mer; Ts have been removed before.
awk -F"\t" '{if($3 == 1){print $1}}' 18.tab | awk -F"separatortoto" '{print $1"\n"$2}' - > Unique-18mers-noT.fasta

##Filter for up to one MM as a 20mer
../QueryInGenomeWithMMinWindow.pl Celegans-mRNA_WB270.fasta Unique-18mers-noT.fasta 20 1 2 19 | awk -F"\t" '{if($3==1){print $1"\n"$2}}' > Unique18mers-1MM-noT.fasta

##Obtain positions of uniquely mapped by bwa
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Unique18mers-1MM-noT.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Unique18mers-1MM-noT.fasta > Cel_exo_1MM_piRNA-seqs.sam
grep "XT:A:U" Cel_exo_1MM_piRNA-seqs.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_exo_1MM_piRNA-seqs.bed

###2 mismatches by blast strategy only
##Create DB and filter for alignements larger than 17
#2 MM

##After installing blast, it is required to create a database of the genome to map against
makeblastdb -in Celegans-mRNA_WB270.fasta -dbtype 'nucl' -title c_elegans_exome.WS270

##Blast filtered unique 18 mer sequences using "blast short algorithm"
blastn -db Celegans-mRNA_WB270.fasta -query Unique18mers-1MM-noT.fasta -task "blastn-short" -out Blast18mers-1MM.tab -outfmt 6 -evalue 1 -num_threads 20

##Get names of seqs that are unique in a set that matches at least 18 bp and start the aligment before the third position
awk -F"\t" '{if((($3*$4)/100)>17){if($7 < 3){array[$1]++}}} END{for (key in array){if(array[key]==1){print key}}}' Blast18mers-1MM.tab > Blast_18mers_filt-2MM.txt

##Parse output to create a file with the indexes of the unique sequences
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' Unique18mers-1MM-noT.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - Blast_18mers_filt-2MM.txt > Cel_exo_2MM_piRNA-seqs.fasta

##Map unique to obtain position
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Cel_exo_2MM_piRNA-seqs.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Cel_exo_2MM_piRNA-seqs.fasta > Cel_exo_2MM_piRNA-seqs.sam
grep "XT:A:U" Cel_exo_2MM_piRNA-seqs.sam| awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_exo_2MM_piRNA-seqs.bed

###3 mismatches by blast strategy only
## Use previous alignement and filter from it
#3 MM

##Get names of seqs that are unique in a set that matches at least 18 bp and start the aligment before the third position
awk -F"\t" '{if((($3*$4)/100)>16){if($7 < 4){array[$1]++}}} END{for (key in array){if(array[key]==1){print key}}}' Blast18mers-1MM.tab > Blast_18mers_filt-3MM.txt

##Parse output to create a file with the indexes of the unique sequences
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' Unique18mers-1MM-noT.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - Blast_18mers_filt-3MM.txt > Cel_exo_3MM_piRNA-seqs.fasta

##Map unique to obtain position
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Cel_exo_3MM_piRNA-seqs.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Cel_exo_3MM_piRNA-seqs.fasta > Cel_exo_3MM_piRNA-seqs.sam
grep "XT:A:U" Cel_exo_3MM_piRNA-seqs.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_exo_3MM_piRNA-seqs.bed

###4 mismatches by 16-mer and blast strategy
## Use previous unique mer fasta and convert it into 16 mer (2MM-16-2MM)
#4 MM
awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-19); i++){print name"_"i"separatortoto"substr($0,i,20)"\n"substr($0,i+2,16)}}}' Cel_exo_3MM_piRNA-seqs.fasta > 3MM-CDS-16mers.fasta

##Look for uniqueness in 18 mers with perl script (hash method)
../QueryInGenomeWithMMinWindow2.pl Celegans-mRNA_WB270.fasta 3MM-CDS-16mers.fasta 16 0 1 16 16.tab

##Convert into fasta files sequences that appear only once, place it as 20mer; Ts have been removed before.
awk -F"\t" '{if($3 == 1){print $1}}' 16.tab | awk -F"separatortoto" '{print $1"\n"$2}' - > Unique-16mers-noT-3MM.fasta

##Blast filtered unique 16 mer sequences using "blast short algorithm"
blastn -db Celegans-mRNA_WB270.fasta -query Unique-16mers-noT-3MM.fasta -task "blastn-short" -out Blast16mers-3MM.tab -outfmt 6 -evalue 1 -num_threads 20

##Get names of seqs that are unique in a set that matches at least 16 bp and start the aligment before the third position
awk -F"\t" '{if((($3*$4)/100)>15){if($7 < 5){array[$1]++}}} END{for (key in array){if(array[key]==1){print key}}}' Blast16mers-3MM.tab > Blast_16mers_filt-4MM.txt

##Parse output to create a file with the indexes of the unique sequences
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' Unique-16mers-noT-3MM.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - Blast_16mers_filt-4MM.txt > Cel_exo_4MM_piRNA-seqs.fasta

##Map unique to obtain position
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Cel_exo_4MM_piRNA-seqs.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Cel_exo_4MM_piRNA-seqs.fasta > Cel_exo_4MM_piRNA-seqs.sam
grep "XT:A:U" Cel_exo_4MM_piRNA-seqs.sam| awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_exo_4MM_piRNA-seqs.bed

###4 mismatches by 15-mer sided to right and then a final blast
## Use previous unique mer fasta and convert it into 15 mer (5MM-15 unique)
#5 MM
awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-19); i++){print name"_"i"separatortoto"substr($0,i,20)"\n"substr($0,i+5,15)}}}' Cel_exo_4MM_piRNA-seqs.fasta > 4MM-CDS-R15mers.fasta

##Look for uniqueness in 15 mers with perl script (hash method)
../QueryInGenomeWithMMinWindow2.pl Celegans-mRNA_WB270.fasta 4MM-CDS-R15mers.fasta 15 0 1 15 15.tab

##Convert into fasta files sequences that appear only once, place it as 20mer; Ts have been removed before.
awk -F"\t" '{if($3 == 1){print $1}}' 15.tab | awk -F"separatortoto" '{print $1"\n"$2}' - > Unique-R15mers-noT-4MM.fasta

##Blast filtered unique 16 mer sequences using "blast short algorithm"
blastn -db Celegans-mRNA_WB270.fasta -query Unique-R15mers-noT-4MM.fasta -task "blastn-short" -out Blast15mers-4MM.tab -outfmt 6 -evalue 1 -num_threads 20

##Get names of seqs that are unique in a set that matches at least 16 bp and start the aligment before the third position
awk -F"\t" '{if((($3*$4)/100)>14){if($7 < 6){array[$1]++}}} END{for (key in array){if(array[key]==1){print key}}}' Blast15mers-4MM.tab > Blast_16mers_filt-5MM.txt

##Parse output to create a file with the indexes of the unique sequences
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' Unique-R15mers-noT-4MM.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - Blast_16mers_filt-5MM.txt > Cel_exo_5MM_piRNA-seqs.fasta

#Align sequences back to the genome and convert map to a bed file
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Cel_exo_5MM_piRNA-seqs.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Cel_exo_5MM_piRNA-seqs.fasta > Cel_exo_5MM_piRNA-seqs.sam
grep "XT:A:U" Cel_exo_5MM_piRNA-seqs.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_exo_5MM_piRNA-seqs.bed



###Please note that any of the fasta files were filtered for GC content; in case that that GC filtering is needed, it can be done in the fasta files with the following command:
#: awk '{if($0 ~ />/){name=$0;}else{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq=$0;k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; GCcontent=(sumC+sumG)/k*100; if((GCcontent >= 30) && (GCcontent <= 45)){print name"_GC-"GCcontent"\n"seq}}}' File.fasta > GCFilt.fasta
## Where 30 is lower boundary and 45 is upper boundary, and File.fasta is input fasta and GCfilt.fasta is name of output fasta.
##Also, this file can be mapped once again with bwa and converted in bed with following command:
#: bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa GCFilt.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - GCFilt.fasta > GCFilt.sam

#: grep "XT:A:U" GCFilt.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > GCFilt.bed




