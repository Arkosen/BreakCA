# BreakCA
Detecting small and medium Indels using ATAC and ChIP-seq reads.

# Introduction: 
Sequencing reads that span variant breakpoints are “chimeric” in local alignment because they either appear to result from the fusion of two sequences or contain insertions/deletions within the read sequence and are discordant from the reference genome. We utilize these chimeric reads in a machine-learning framework to detect indels within ATAC and ChIP-seq peaks

# Methods: 
Step wise implementation of BreakCA is as follows. All samples are aligned using BWA-MEM using default params. 

1. Get reads overlapping peaks.
Rscript --vanilla  ~/BreakCA/bin/get_reads_from_bam.R input.bam peaks.bed reads.tsv

2. Get read pileups; allow for correction for mis-alignments
samtools mpileup -B -f genome.fa -l peaks.bed input.bam | gzip > read.pileups.gz

3. Convert pileups to readable files
less read.pileups.gz| awk -v OFS='\t' '{print $1,$2,$5}' > read.pileup

4. Get insertion positions
less read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /\+[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > insertion.pileups

5. Get deletion containing positions
less read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /-[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > deletion.pileups

6. Counts reads at each bp position within peaks
Rscript --vanilla ~/BreakCA/bin/count_reads_per_base.R reads.tsv insertion.pileups deletion.pileups read.pileup sc.tsv counts.tsv

7. Add clipping information including information content and clip length
Rscript --vanilla ~/BreakCA/bin/get_clipping_information.R sc.tsv sc_w_seq.tsv clip.info.txt

8. Calculate posterior mean and standard deviation for clipped positions and output contig assembly regions as bed files.
Rscript --vanilla ~/BreakCA/bin/calculate_posterior.R counts.tsv posteriors.tsv contig_assembly.bed clip.info.txt all.positions.tsv

9. Run contig assembly

10. get read ID
samtools view -L contig_assembly.bed input.bam | cut -f1 > read.Ids.txt

11. If pe (pe=1) get fastq for pair-end format or get single end (pe=0) format
# for first pair (implemented in Linux)
if (( $pe > 0 ))
then seqtk subseq forward_reads.fa read.Ids.txt  > r1.fastq
else seqtk subseq forward_reads.fa read.Ids.txt > unpaired.fastq
fi
# for second pair (implemented in linux)
if (( $pe > 0 ))
then seqtk subseq reverse_reads.fa read.Ids.txt  > r2.fastq
fi

12. Run SPAdes contig assembly 
# implemented in linux
if (( $pe > 0 ))
then spades.py -t 12 -k 21,33 --only-assembler -o spades -1 r1.fastq -2 r2.fastq
else spades.py -t 12 -k 21,33 --only-assembler -o spades -s unpaired.fastq
fi

13. Align SPAdes
# implemented in linux
perl ~/BreakCA/bin/fasta_to_fastq.pl spades/contigs.fasta > contigs.fq
bwa mem -t 12 genome.fa contigs.fq | samtools view -S -b -h -F 4 - > contigs.bam
samtools sort contigs.bam contigs.sorted
samtools index contigs.sorted.bam

14. Find contig supported base positions
results=posteriors.tsv
contigs=contigs.sorted.bam
Rscript --vanilla ~/BreakCA/bin/contig_support_wrapper.R $contigs $results contig.supp.txt

15. Predefine regions to test, 20bps non-overlapping windows
Rscript --vanilla ~/BreakCA/bin/predefine_windows.R peaks.bed windows.bed

16. Prepare contigs for machine-learning
Rscript --vanilla ~/BreakCA/bin/prepare_dataset.R all.positions.tsv windows.bed contig.supp.txt classifier_id_frame.csv classifier_input.tsv

17. Make prediction using logistic regression, models can be created using build_randomForest.R script under ~/BreakCA/misc, we also provide pe and se models created using GM12878 ATAC-seq and ChIP-seq reads.
Rscript --vanilla ~/BreakCA/bin/make_predictions.R classifier_input.tsv $model prediction.txt

# The feature map building scripts can be run using breakCA.bash shell script in linux
Usage= ./breakCA.bash -a <path to R> -b <.bam> -p <.bed>  -f <read1.fq.gz> -r <read2.fq.gz> -o <output directory> -s <1=paired end, 0=single end > -g <fasta file for genome, directory should contain bwa index> -t <# of threads>

# Prediction on feature map can be can be run using predict.bash shell script in linux
Usage= ./predict.bash -a <path to R> -o <output directory> -w < peaks > -m < model >

# Models are available for download in ~/BreakCA/models

