# BreakCA
Detecting small and medium Indels using ATAC and ChIP-seq reads.

# Introduction: 
Sequencing reads that span variant breakpoints are “chimeric” in local alignment because they either appear to result from the fusion of two sequences or contain insertions/deletions within the read sequence and are discordant from the reference genome. We utilize these reads in a machine-learning framework to detect indels within ATAC and ChIP-seq peaks

# Methods: 
Step wise implementation of BreakCA is as follows. All samples are aligned using BWA-MEM using default params. 

# 1. Get reads overlapping peaks.
Rscript --vanilla  ~/BreakCA/bin/get_reads_from_bam.R input.bam peaks.bed reads.tsv

# 2. Get read pileups; allow for correction for mis-alignments
samtools mpileup -B -f genome.fa -l peaks.bed input.bam | gzip > read.pileups.gz

# 3. Convert pileups to readable files
less read.pileups.gz| awk -v OFS='\t' '{print $1,$2,$5}' > read.pileup

# 4. Get insertion positions
less read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /\+[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > insertion.pileups

# 5. Get deletion containing positions
less read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /-[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > deletion.pileups

# 6. Counts reads at each bp position within peaks
Rscript --vanilla ~/BreakCA/bin/count_reads_per_base.R reads.tsv insertion.pileups deletion.pileups read.pileup sc.tsv counts.tsv

# 7. Add clipping information including information content and clip length
Rscript --vanilla ~/BreakCA/bin/get_clipping_information.R sc.tsv sc_w_seq.tsv clip.info.txt

# 8. Calculate posterior mean and standard deviation.
Rscript --vanilla ~/BreakCA/bin/calculate_posterior.R counts.tsv posteriors.tsv clip.info.txt all.positions.tsv

# 9. Predefine regions to test, 20bps non-overlapping windows
Rscript --vanilla ~/BreakCA/bin/predefine_windows.R peaks.bed windows.bed

# 10. Prepare contigs for machine-learning
Rscript --vanilla ~/BreakCA/bin/prepare_dataset.R all.positions.tsv windows.bed classifier_id_frame.csv classifier_input.tsv

# 11. Make prediction using logistic regression, 
Models can be created using build_randomForest.R script under ~/BreakCA/misc. We also have pe and se models created using GM12878 ATAC-seq and ChIP-seq reads available on request.

Rscript --vanilla ~/BreakCA/bin/make_predictions.R classifier_input.tsv model.rda prediction.txt

# The feature map building scripts can be run using breakCA.bash shell script in linux
Usage= ./breakCA.bash -a <path to R> -b <.bam> -p <.bed> -o <output directory> -g <fasta file for genome> 

# Prediction on feature map can be can be run using predict.bash shell script in linux
Usage= ./predict.bash -a <path to R> -o <output directory> -w < peaks > -m < model >


