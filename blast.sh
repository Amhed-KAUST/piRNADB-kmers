##After installing blast, it is required to create a database of the genome to map against
makeblastdb -in c_elegans.PRJNA13758.WS271.genomic_softmasked.fa -dbtype 'nucl' -title c_elegans.PRJNA13758.WS271

##Blast filtered unique sequences using "blast short algorithm"
blastn -db ../db/c_elegans.PRJNA13758.WS270.genomic.fa -query Unique14mers-20mers.3MM-3-17.fasta -task "blastn-short" -out Unique14mers-20mers.3MM-3-17.tab -outfmt 6 -evalue 1 -num_threads 20

##Parse output to create a file with the indexes of the unique sequences
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' Unique14mers-20mers.3MM-3-17.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - indexes.names > Unique14mers-20mers.3MM-3-17.blast15aF.fasta

#Align sequences back to the genome and convert map to a bed file
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Unique14mers-20mers.3MM-3-17.blast15aF.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Unique14mers-20mers.3MM-3-17.blast15aF.fasta > 3MM-15aln.sam
