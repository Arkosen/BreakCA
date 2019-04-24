# BreakCA
Detecting recurrent mutated windows using ATAC and ChIP-seq reads.

# Introduction: 
Sequencing reads that span variant breakpoints are “chimeric” in local alignment because they either appear to result from the fusion of two sequences or contain insertions/deletions within the read sequence and are discordant from the reference genome. We utilize these reads in a scalable machine-learning framework to detect indels within ATAC and ChIP-seq peaks.

# System requirements
R(version= 3.3.1), samtools (Version: 0.1.19-44428cd), GATK (version= 3.7-0-gcfedb67) or GATK (version= 4.1.1.0)

Software was only tested on above versions. No non-standard hardware required

# Installation guide
#CRAN packages

install.packages(c('data.table', 'plyr', 'dplyr', 'pbapply', 'readr', 'reshape', 'rmutils', 'lattice', 'stringr', 'mlr'))

#Bioconductor

source("https://bioconductor.org/biocLite.R");
biocLite(c("biovizBase", "rtracklayer", "Rsamtools", "BSgenome.Hsapiens.UCSC.hg19", "GenomicAlignments", "VariantAnnotation")

# Instruction to run on data
All samples are aligned using BWA-MEM using default params. All read with MAPQ>30 is kept for analysis.

1. Get reads overlapping peaks:

Rscript --vanilla  ~/BreakCA/bin/get_reads_from_bam.R input.bam peaks.bed reads.tsv

2. Get read pileups; allow for correction for mis-alignments:

samtools mpileup -B -f genome.fa -l peaks.bed input.bam | gzip > read.pileups.gz

3. Convert pileups to readable files:

less read.pileups.gz| awk -v OFS='\t' '{print $1,$2,$5}' > read.pileup

4. Get insertion positions:

less read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /\+[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > insertion.pileups

5. Get deletion containing positions:

less read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /-[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > deletion.pileups

6. Counts reads at each bp position within peaks:

Rscript --vanilla ~/BreakCA/bin/count_reads_per_base.R reads.tsv insertion.pileups deletion.pileups read.pileup sc.tsv counts.tsv

7. Add clipping information including information content and clip length:

Rscript --vanilla ~/BreakCA/bin/get_clipping_information.R sc.tsv sc_w_seq.tsv clip.info.txt

8. Calculate posterior mean and standard deviation:

Rscript --vanilla ~/BreakCA/bin/calculate_posterior.R counts.tsv posteriors.tsv clip.info.txt all.positions.tsv

9. Predefine regions to test, 20bps non-overlapping windows (implemented for hg19):

Rscript --vanilla ~/BreakCA/bin/predefine_windows.R peaks.bed windows.bed

10. Prepare for predictions:

Rscript --vanilla ~/BreakCA/bin/prepare_dataset.R all.positions.tsv windows.bed classifier_id_frame.csv classifier_input.tsv

optional (add QD i.e. QualDepth from GATK):

Rscript --vanilla ~/BreakCA/bin/add_QD.R classifier_input.tsv gatk.indels.vcf windows.bed classifier_input_w_QD.tsv 

11. Make prediction:

Models can be created using build_randomForest.R script under ~/BreakCA/misc. This is a time-consuming step but needs to run only once

Rscript --vanilla ~/BreakCA/misc/build_randomForest.R training.txt model.rda

Rscript --vanilla ~/BreakCA/bin/make_predictions.R classifier_input.tsv model.rda prediction.txt 

The feature map building scripts can be run using breakCA.bash shell script in linux:

Usage= ./breakCA.bash -a path to R -b .bam -p .bed -o output directory -g fasta file for genome

Prediction on feature map can be can be run using predict.bash shell script in linux:

Usage= ./predict.bash -a path to R -o output directory -w peaks -m model -g 0

Note: when g=0 QD from GATK is not added, when g=1 QD is added provided the output directory contain file named gark.indels.vcf

Note: We provide constructed training and testing feature map used for GM12878 cell lines as a zip file.

# Expected output 
text file named "prediction.txt" with rownames= "chr_start_end" (e.g. chr1_100_120) and 3 columns c("prob.1KG" "prob.N" "response"). "prob_1KG"= probability that window contain indel which look like a 1000Genomes Indel detected GM12878 ATAC/ChIP-seq dataset (training), "prob.N"= probability its not an indel containing windows, "response"= 1KG or N
