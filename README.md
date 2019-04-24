# BreakCA
A method to detect indels using ATAC-seq and ChIP-seq reads.

# Introduction
BreakCA is short for “Breaks in Chromatin Accessible” regions. It is a software package that detects indels in small (20bp) windows using ChIP-seq or ATAC-seq reads. BreakCA exploits properties of mapped reads such as gaps in alignments (i.e. insertions or deletions) and clipping at read ends.  BreakCA collects 16 features from mapped reads and uses a random forest to identify windows that contain indels. BreakCA is described in the following preprint: https://www.biorxiv.org/content/10.1101/605642v1.abstract

# System requirements
BreakCA requires only a standard computer with enough RAM to support the operations performed by a user. We recommend a computer with the following specs:

* RAM: 12+ GB
* CPU: 4+ cores

# Installation guide
BreakCA is a collection of Rscripts and doesnot require extra-installation steps. Just download the scripts to your home directory and run. To run it requires the following packages installed on R. BreakCA is designed to be run in a unix/linux shell.

## Software versions

BreakCA was tested using the following software versions:
 * R (version= 3.3.1)
 * samtools (version: 0.1.19-44428cd)
 * GATK (version= 3.7-0-gcfedb67) or GATK (version= 4.1.1.0) (OPTIONAL)

## CRAN packages

```R
install.packages(c('data.table', 'plyr', 'dplyr', 'pbapply', 'readr', 'reshape', 'rmutil', 'lattice', 'stringr', 'mlr'))
```

## Bioconductor

```R
source("https://bioconductor.org/biocLite.R");
biocLite(c("biovizBase", "rtracklayer", "Rsamtools", "BSgenome.Hsapiens.UCSC.hg19", "GenomicAlignments", "VariantAnnotation")
```

## Installation time
The installation time depends on the installation time for the required packages, but should only take a matter of minutes.

## Package versions
BreakCA was tested with the following package versions:
* data.table	1.10.0,
* plyr	1.8.4,
* dplyr	0.7.4,
* pbapply	1.3.3,
* readr	1.1.1,
* reshape	0.8.6,
* rmutil	1.1.0,
* lattice	0.20.33,
* stringr	1.3.1,
* mlr	2.12.1,
* biovizBase	1.22.0,
* rtracklayer	1.34.2,
* Rsamtools	1.26.1,
* BSgenome.Hsapiens.UCSC.hg19	1.4.0,
* GenomicAlignments	1.10.0,
* VariantAnnotation	1.20.2

# Running BreakCA

BreakCA takes a BAM file as input and a set of ATAC-seq and ChIP-seq peaks in BED format. 

The BAM file should contain reads that are aligned using BWA-MEM with default paramaters. All reads with MAPQ>30 are kept for analysis.

BreakCA can be run using a pair of wrapper scripts or by running a series of commands corresponding to each step of the analysis pipeline.

## Using wrapper scripts

We provide a pair of bash scripts to run BreakCA. 

The `breakCA.bash` script builds the feature table that is used to make predictions.
```bash
./breakCA.bash -a path_to_R -b path_to_aligned_reads.bam -p path_to_peaks.bed -o output_directory -g fasta_file_for_genome.fa
```

Predictions can be made from the feature table using the `predict.bash` shell script:
```bash
./predict.bash -a path_to_R -o output_directory -w path_to_peaks.bed -m model -g 0
```

The `-g` flag indicates whether QD from GATK should be used in predictions. When `-g 0` QD from GATK is not added. When `-g 1` QD is added, provided the output directory contains a file named `gatk.indels.vcf` which can be generated for each sample using the GATK HaplotypeCaller.

## Using individual commands

BreakCA can be run as as series of separate shell commands and scripts as follows. 

1. Get reads overlapping peaks:
   ```bash
   Rscript --vanilla  BreakCA/bin/get_reads_from_bam.R input.bam peaks.bed reads.tsv
   ```

2. Get read pileups; allow for correction for mis-alignments:
   ```bash
   samtools mpileup -B -f genome.fa -l peaks.bed input.bam | gzip > read.pileups.gz
   ```

3. Convert pileups to more-easily readable files:
    ```bash
    gunzip -c read.pileups.gz | awk -v OFS='\t' '{print $1,$2,$5}' > read.pileup
    ```

4. Get insertion positions:
    ```bash
    gunzip -c read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /\+[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > insertion.pileups
    ```

5. Get deletion containing positions:
    ```bash
    gunzip -c read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /-[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > deletion.pileups
    ```

6. Counts reads at each base position within peaks:
    ```bash
    Rscript --vanilla BreakCA/bin/count_reads_per_base.R reads.tsv insertion.pileups deletion.pileups read.pileup sc.tsv counts.tsv
    ```

7. Add clipping information including information content and clip length:
    ```bash
    Rscript --vanilla BreakCA/bin/get_clipping_information.R sc.tsv sc_w_seq.tsv clip.info.txt
    ```

8. Calculate posterior means and standard deviations:
    ```bash
    Rscript --vanilla BreakCA/bin/calculate_posterior.R counts.tsv posteriors.tsv clip.info.txt all.positions.tsv
    ```

9. Predefine regions to test. These are 20bp non-overlapping windows (implemented for hg19):
    ```bash
    Rscript --vanilla BreakCA/bin/predefine_windows.R peaks.bed windows.bed
    ```

10. Prepare data table containing features for predictions:
    ```bash
    Rscript --vanilla BreakCA/bin/prepare_dataset.R all.positions.tsv windows.bed classifier_id_frame.csv classifier_input.tsv
    ```

11. Add QD (i.e. QualDepth) from GATK (OPTIONAL):
    ```bash
    Rscript --vanilla BreakCA/bin/add_QD.R classifier_input.tsv gatk.indels.vcf windows.bed classifier_input_w_QD.tsv 
    ```
    
12. Train the Random Forest classifier. Note this is a time-consuming step but it only needs to be run once:
    ```bash
    Rscript --vanilla BreakCA/misc/build_randomForest.R training.txt model.rda
    ```

13. Predict indels using the trained Random Forest:
    ```bash
    Rscript --vanilla BreakCA/bin/make_predictions.R classifier_input.tsv model.rda prediction.txt 
    ```
 
Note: We provide feature tables from the GM12878 cell line for training and testing as a zip file.

# Expected output 

The `make_predictions.R` script outputs a text file `prediction.txt` with rownames = "chr_start_end" (e.g. chr1_100_120) and 3 columns ("prob.1KG", "prob.N", "response"). "prob_1KG" is BreakCA score, and can be roughly interpreted as the probability that a window contains an indel resembling those found in GM12878 ATAC/ChIP-seq training dataset. "prob.N" is the probability that the window does not contain an indel, "response"= 1KG or N
