#!/bin/sh
# help
if [ "$1" == "-h" ] ; then
    echo "Usage= ./breakCA.bash -a <path to R> -b <.bam> -p <.bed>  -f <read1.fq.gz> -r <read2.fq.gz> -o <output directory> -s <1=paired end, 0=single end > -g <fasta file for genome, directory should contain bwa index> -t <# of threads>"
    exit 0
fi

# running with options
while getopts ":a:b:p:f:r:o:s:g:t:" opt
   do
     case $opt in
		a ) path=$OPTARG;;
        b ) bam=$OPTARG;;
        p ) peaks=$OPTARG;;
        f ) forward_reads=$OPTARG;;
        r ) reverse_reads=$OPTARG;;
		o ) output=$OPTARG;;
		s ) pe=$OPTARG;;
		g ) genome=$OPTARG;;
		t ) threads=$OPTARG;;
     esac
done

# make output directory
mkdir -p $output

# get reads as a tsv file
$path/Rscript --vanilla  ~/BreakCA/bin/get_reads_from_bam.R $bam $peaks $output/read.tsv

# get read pileups; allow for correction for mis-alignments
samtools mpileup -B -f $genome -l $peaks $bam | gzip > $output/read.pileups.gz

# get bases written into a file to use during counting
less $output/read.pileups.gz| awk -v OFS='\t' '{print $1,$2,$5}' > $output/read.pileup

# get insertion containing positions
less $output/read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /\+[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > $output/insertion.pileups

# get deletion containing positions
less $output/read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /-[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > $output/deletion.pileups

# counts reads at each bp position within peaks
$path/Rscript --vanilla ~/BreakCA/bin/count_reads_per_base.R $output/read.tsv $output/insertion.pileups $output/deletion.pileups $output/read.pileup $output/sc.tsv $output/counts.tsv

  # add clipping information including information content and clip length
$path/Rscript --vanilla ~/BreakCA/bin/get_clipping_information.R $output/sc.tsv $output/sc_w_seq.tsv $output/clip.info.txt

  # calculate posterior mean and standard deviation for clipped positions
$path/Rscript --vanilla ~/BreakCA/bin/calculate_posterior.R $output/counts.tsv $output/posteriors.tsv $output/contig_assembly.bed $output/clip.info.txt $output/all.positions.tsv

  # run contig assembly
 echo "running contig assembly"

 bed=$output/contig_assembly.bed

  #get read name
samtools view -L $bed $bam | cut -f1 > $output/read.Ids.txt

# if pe>0 get fastq for pair-end reads or get single end
# for first pair
if (( $pe > 0 ))
then seqtk subseq $forward_reads $output/read.Ids.txt  > $output/r1.fastq
else seqtk subseq $forward_reads $output/read.Ids.txt > $output/unpaired.fastq
fi

# for second pair
if (( $pe > 0 ))
then seqtk subseq $reverse_reads $output/read.Ids.txt  > $output/r2.fastq
fi

 # run SPAdes
if (( $pe > 0 ))
then spades.py -t $threads -k 21,33 --only-assembler -o $output/spades -1 $output/r1.fastq -2 $output/r2.fastq
else spades.py -t $threads -k 21,33 --only-assembler -o $output/spades -s $output/unpaired.fastq
fi

 # align SPAdes
perl ~/BreakCA/bin/fasta_to_fastq.pl $output/spades/contigs.fasta > $output/contigs.fq
bwa mem -t $threads $genome $output/contigs.fq | samtools view -S -b -h -F 4 - > $output/contigs.bam
samtools sort $output/contigs.bam $output/contigs.sorted
samtools index $output/contigs.sorted.bam

# find contig supported base positions
results=$output/posteriors.tsv
contigs=$output/contigs.sorted.bam
$path/Rscript --vanilla ~/BreakCA/bin/contig_support_wrapper.R $contigs $results $output/contig.supp.txt

# run varscan2:: only required for precision recall testing
#samtools view -b -h -L $peaks $bam > $output/peaks.bam
#samtools mpileup -B -f $genome $output/peaks.bam | varscan mpileup2indel --min-coverage 10 --p-value 1 --output-vcf 1 > $output/indels.vcf

echo "done"
