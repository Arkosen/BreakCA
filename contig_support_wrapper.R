args <- commandArgs(trailingOnly = TRUE)
contigs<- args[1]
results<- args[2]
filename.contigs<- args[3]
print(args)

# load library
library(GenomicAlignments)
library(rtracklayer)
library(Rsamtools)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg19)
library(lattice)
library(plyr)
library(dplyr)
library(data.table)
library(reshape)
library(rmutil)
library(lattice)
library(stringr)
library(readr)

# read contigs bam
what = c("seq","qual")
param = ScanBamParam(what=what)
bam = readGAlignments(contigs, param = param)
dat=as.data.frame(bam)
print("read contig file")

#get features
feature= c("H","S","I","D")
print(paste("features types are",feature, sep=" "))

# grep feature containing contigs
grep.func= function(i){
   df= dat[grep(feature[i], dat$cigar),]
   return(df)
}
df= do.call(rbind, lapply(1:length(feature),grep.func))
df= df[!duplicated(df),]
df$id = paste(df$seqnames, df$start, df$end, df$cigar,sep = ";")

# get contig cigar
cigar= df$cigar
ref= cigarRangesAlongReferenceSpace(cigar = cigar)
query= cigarRangesAlongQuerySpace(cigar=cigar)
opt= explodeCigarOps(cigar, ops=CIGAR_OPS)
print("CIGAR ranges and form collected")

# get ranges along the reference space in cigar
get.ranges= function(i){
   # get feature information
   r= as.data.frame(ref[[i]])
   q= as.data.frame(query[[i]])
   r$width.query.space= q$width
   d= df[i,]
   o= opt[[i]]
   
   #parse feature information
   myfunc= function(n){
      seq= as.character(d$seq)   
      pos= q[n,]
      sub.seq= substr(seq, start=pos$start, stop=pos$end)
      pos$seq= sub.seq
      pos$seqlen= nchar(sub.seq)
      return(pos)
   }
   q= do.call(rbind, lapply(1:nrow(q), myfunc))
   r$seq.query.space=q$seq
   r$start= r$start + (d$start)
   r$end= r$end + (d$start)
   r$opt= o
   r$id= paste(d$seqnames, d$start, d$end, d$cigar,sep = ";")
   r$cigar=d$cigar
   r$space= d$seqnames
   r
}
cigar.list= llply(1:length(cigar), get.ranges,.progress = progress_text(char = "=") )
cigar.coord= as.data.frame(rbindlist(cigar.list))
print("parsed coordinates in reference space")

# find deletion contigs
deletion<- subset(cigar.coord, cigar.coord$opt== "D")
print("found deletion containing contigs")

# find insertion, soft-clipped or hard-clipped contigs and correcting the reference position
type= c("S","I", "H")
correct= function(i){
   dl= cigar.coord[cigar.coord$opt %in% type[i],]
   dl$end= dl$start
   dl$start= dl$start-1
   dl$end= dl$end-1
   return(dl)
}
other= do.call(rbind, lapply(1:length(type), correct))
print("corrected reference position")

# combine deletion and other contigs
mat= rbind.fill(deletion,other)
print("contig features parsed in reference space")

# overlap with contigs
gr= as(as(mat, "RangedData"), "GRanges")

# load results
results<- as.data.frame(fread(results))
results$seqnames<- sapply(strsplit(results$bpos, ":"), "[", 1)
results$start<- as.numeric(sapply(strsplit(results$bpos, ":"), "[", 2))
results$end<- results$start
query= makeGRangesFromDataFrame(results, keep.extra.columns = T)

# find nearest
dist<- distanceToNearest(query, gr, ignore.strand=FALSE)
cs<- as.data.frame(query[queryHits(dist)])
con<- as.data.frame(gr[subjectHits(dist)])
con$bpos<- cs$bpos
print("getting closest contig features near break position")

#match and transfer information
con$dist<- mcols(dist)$distance
final<- con[con$dist <=10,]
final$seq<- df$seq[match(final$id, df$id)]
write.table(final,filename.contigs, row.names=F, sep="\t")
print("finding contig features within 10bps of break positions")
