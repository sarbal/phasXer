---
title: "Phasing the X"
author: "Sara Ballouz"
date: "November 7, 2019"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


# Aligning single-cell 
## Using STARsolo 
Example with carrier 1 data. Repeat this for all samples. Change the read names (in your script).  
```{}
STAR \
--runThreadN 15 \
--soloType Droplet \
--genomeDir /data/genomes/GRCh38_Gencode31/ \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--soloCBwhitelist /data/genomes/cellranger_barcodes/3M-february-2018.txt \
--readFilesIn SCGC-GILL-JG-03_S1_L001_R2_001.fastq.gz,SCGC-GILL-JG-03_S1_L002_R2_001.fastq.gz,SCGC-GILL-JG-03_S1_L003_R2_001.fastq.gz,SCGC-GILL-JG-03_S1_L004_R2_001.fastq.gz SCGC-GILL-JG-03_S1_L001_R1_001.fastq.gz,SCGC-GILL-JG-03_S1_L002_R1_001.fastq.gz,SCGC-GILL-JG-03_S1_L003_R1_001.fastq.gz,SCGC-GILL-JG-03_S1_L004_R1_001.fastq.gz \
--readFilesCommand zcat \
--soloCellFilter  CellRanger2.2 8000 0.99 10 \
--soloFeatures Gene Velocyto \
--soloBarcodeReadLength 28 \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
```

## Using Cellranger
This is the script that the single-cell core (Jon Preall) uses to run cellranger on their server. 
Example for control 1, monocytes is below. A clunkier version of this is at the bottom, which can be run on our own servers. 
```{}
#!/bin/sh
#$ -cwd
#$ -v PATH,LD_LIBRARY_PATH

BASEDIR=$(pwd)
#mkdir -p $BASEDIR/mkfastq || exit 1
mkdir -p $BASEDIR/count || exit 1
#mkdir -p $BASEDIR/aggr || exit 1

#Stage 1 - make fastq files


#run this script from within the desired experiment home directory
#script will create a mkfastq directory in which to deposit the output of the mkfastq pipline
#Before running, make sure that a valid SampleSheet.csv file is located inside of the home directory

#Specify RUNID
RUNID=190823_NB501555_0625_AH2HYJBGXC
FLOWCELL=$(echo $RUNID | grep -o '.\{9\}$')

#path to basecalls
BCL=/seq/Illumina_runs/NextSeqData/NextSeqOutput/$RUNID

#path to basecalls
BCL=/seq/Illumina_runs/NextSeqData/NextSeqOutput/$RUNID

#Name your experiment:
EXPNAME=Gillis_05

#cd $BASEDIR/mkfastq
#cellranger mkfastq \
#       --run=$BCL \
#       --samplesheet=$BASEDIR/SampleSheet.csv \
#       --jobmode=sge \
#        || exit 1
#cd $RUNDIR
#mail -s "mkfastq pipeline completed for $RUNID" "$USER@cshl.edu" || exit 1

#Select your transcriptome by decommenting the appropriate line:
#TRANSCRIPTOME=/seq/CellRanger/references/Diermeier
#TRANSCRIPTOME=/seq/CellRanger/references/IMGT_Mouse_Reference_for_VDJ
#TRANSCRIPTOME=/seq/CellRanger/references/MAPSeq_files
#TRANSCRIPTOME=/seq/CellRanger/references/mm10-2.1.0.premrna
#TRANSCRIPTOME=/seq/CellRanger/references/mm10_MAPSeq
#TRANSCRIPTOME=/seq/CellRanger/references/mouse_vdj_reference_package_Ensembl
#TRANSCRIPTOME=/seq/CellRanger/references/refdata-cellranger-ercc92-1.2.0
TRANSCRIPTOME=/seq/CellRanger/references/refdata-cellranger-GRCh38-1.2.0
#TRANSCRIPTOME=/seq/CellRanger/references/refdata-cellranger-hg19_and_mm10-1.2.0
#TRANSCRIPTOME=/seq/CellRanger/references/refdata-cellranger-mm10-1.2.0
#TRANSCRIPTOME=/seq/CellRanger/references/refdata-cellranger-mm10-2.1.0
#TRANSCRIPTOME=/seq/CellRanger/references/refdata-cellranger-mm10-3.0.0
#TRANSCRIPTOME=/seq/CellRanger/references/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0
#TRANSCRIPTOME=/seq/CellRanger/references/Zea_Mays

##YOU SHOULDN'T HAVE TO MODIFY ANYTHING BENEATH THIS LINE

##Stage 2: Run Cellranger Count
cd $BASEDIR/count

#Stage 2a - count sample 1
cellranger count --id=Gillis_05_XCGD10xGEX_PBMC \
        --jobmode=sge \
        --transcriptome=$TRANSCRIPTOME \
        --fastqs=$BCL/$FLOWCELL/outs/fastq_path \
        --sample=Gillis_05_XCGD10xGEX_PBMC \
        --project=Gillis_05

echo "CellRanger count complete for $BASEDIR Gillis_05_XCGD10xGEX_PBMC" | mail -s "CellRanger Count Output" "$USER@cs
hl.edu"

```

# Calling variants from single-cell data 
## Using GATK v3 
```{}
picard=/sonas-hs/gillis/hpc/home/sballouz/GenomeAnalysisTK/picard.jar
GATK=/sonas-hs/gillis/hpc/home/sballouz/GenomeAnalysisTK/GenomeAnalysisTK.jar
igvtools=/sonas-hs/gillis/hpc/home/sballouz/IGVTools/igvtools.jar
samtools=/opt/hpc/bin/samtools


### STAR genome 
genome_fa=/sonas-hs/gillis/hpc/home/sballouz/STAR/GRCh38_Gencode25/GRCh38.p7.genome.fa
X_fa=/sonas-hs/gillis/hpc/home/sballouz/STAR/GRCh38_Gencode25/chrX.fa
bam_file=Aligned.sortedByCoord.out.bam


### Cellranger 
genome_fa=/sonas-hs/gillis/hpc/home/sballouz/lyon/XSKEW/sc/10x/v3_chemistry/10k_pbmcs_healthy_donor/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa 
X_fa=/sonas-hs/gillis/hpc/home/sballouz/lyon/XSKEW/sc/10x/v3_chemistry/10k_pbmcs_healthy_donor/refdata-cellranger-GRCh38-3.0.0/fasta/X.fa
bam_file=possorted_genome_bam.bam



samtools view -H ../possorted_genome_bam_sorted_chrX.split.filtered.bam  > header.sam


```


### From STARsolo output 
```{}
chr=chrX

if [ ! -e "$inDir/Aligned.sortedByCoord.out.bam" ]
then
	$samtools sort -@ 10 $inDir/Aligned.out.bam Aligned.sortedByCoord.out
fi

if [ ! -e "$inDir/Aligned.sortedByCoord.out.bai" ]
then
	$samtools index $inDir/Aligned.sortedByCoord.out.bam
fi

```



### From Cellranger output 
```{}

samtools view -h $bam_file.split.filtered.bam  |  grep -E '(^@|CB\:Z\:)' |  samtools view -S -b -o  $bam_file.split.filtered.barcoded.bam  -

samtools index $bam_file.split.filtered.barcoded.bam

```


## Variant calling GATK v3
```{}
$samtools view -b $inDir/Aligned.sortedByCoord.out.bam $chr > $outDir/$chr.bam
$samtools view -b -q 10 $outDir/$chr.bam > $outDir/$chr.filt.bam


echo "Adding read groups to bam file"
java -jar $picardAddOrReplaceReadGroups \
    I=possorted_genome_bamX.filt.bam \
    O=possorted_genome_bamX.rg.bam \
    SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

echo "Marking duplicates"
java -jar $picard MarkDuplicates \
    I=possorted_genome_bamX.rg.bam  \
    O=possorted_genome_bamX.dedupped.bam \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

echo "Splitting and trimming"
java -jar $GATK \
   -T SplitNCigarReads         \
   -R $genome_fa \
   -I possorted_genome_bamX.dedupped.bam \
   -o possorted_genome_bamX.split.filtered.bam \
   -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS


echo "Haplotype calling (?)"
java -jar $GATK \
   -T HaplotypeCaller \
   -R $genome_fa \
   -L X \
   -I possorted_genome_bamX.split.filtered.bam \
   -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 \
   -o possorted_genome_bamX.split.filtered.vcf

echo "Counting SNPs"
java -jar $igvtools count -z 0 -w 1 --bases --strands read   \
   possorted_genome_bamX.split.filtered.bam   \
   possorted_genome_bamX.split.filtered.wig   \
   $X_fa

 

	grep ^X possorted_genome_bamX.split.filtered.vcf  | cut -f2-5 > A
	cut -f1 A > B
	grep -Fw -f B /sonas-hs/gillis/hpc/home/sballouz/lyon/ncbi/dbGAP-11987/C > D
	grep -Fw -f D possorted_genome_bamX.split.filtered.wig  > possorted_genome_bamX.intersect.wig 
	grep -Fw -f D A > possorted_genome_bamX.cut
	rm A B D 

```

## Variant calling GATK v4
```{}
# Need to run STAR with these parameters: 
#run STAR
$STAR  --runThreadN 4 \
 --genomeDir  $genomedir  \
 --readFilesIn $reads \
 --readFilesCommand zcat \
 --outSAMtype BAM SortedByCoord  \
 --quantMode GeneCounts  \
 --twopassMode Basic  \
 --twopass1readsN -1 \
 --outSAMmapqUnique 60
 
 
echo "Adding read groups to bam file"
java -jar $picard AddOrReplaceReadGroups \
    I=$chr.filt.bam \
    O=$chr.rg.bam \
    SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

echo "Marking duplicates"
java -jar $picard MarkDuplicates \
    I=$chr.rg.bam  \
    O=$chr.dedupped.bam \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

echo "Splitting and trimming"
$GATK SplitNCigarReads \
   -R $genome_fa \
   -I chrX.dedupped.bam \
   -O chrX.split.filtered.bam

echo "Haplotype calling (?)"
$GATK HaplotypeCaller \
   -R $genome_fa \
   -L $chr \
   -I  $chr.split.filtered.bam \
   -O  $chr.filtered.vcf \
   --dont-use-soft-clipped-bases \
   --standard-min-confidence-threshold-for-calling 20.0



### these two commands are for the "bulk" version
$GATK ASEReadCounter \
   -R $genome_fa \
   -O $chr.counts.csv \
   -I $chr.split.filtered.bam \
   -V $chr.filtered.vcf

echo "Counting SNPs"
java -jar $igvtools count -z 0 -w 1 --bases --strands read   \
   $chr.split.filtered.bam   \
   $chr.split.filtered.wig   \
   $genome_fa

```

# Getting SNP counts per cell
This step is somewhat inefficient and could be adjusted. 
Idea here is to split the BAM file into single bam files per cell. As all the filtering has been done prior, this step is mostly to get cell specific information. 
```{} 
# Run GATK pipeline 
bash run_skew.sh file

samtools view -H ../possorted_genome_bam_sorted_chrX.split.filtered.bam  > header.sam

# include ID of folder/see script for more details
# Slit bam file into barcode specific samfiles
perl split_bam_barcode.pl 298818

# file <- list of IDS 

# make the wig files from the sam files 
bash make_wig.sh 

# intersect the wig files with known SNPs from dbGAP 
bash intersect_wig.sh


```

# Parsing and generating cell by SNP matrices
```{}
# Set up environment
source("~/useful4.r")

# Load annotations that match the annotations from the mapping
load("~/sballouz/genome/gene_annotations_v29.Rdata")

# These are the PAR regions based on GRCh38 coordinates. 
par1 = 10000:2781479
par2 = 155701382:156030895

# These are the WIG files generated from the BAM file. See the split bam section above. 
input_files = "/mnt/grid/gillis/hpc/data/data/sballouz/human/RNAseq/XCGD/pilot/outs/split/wig_files"

# This the file with the variant positions (VCF) called from the RNA-seq data. See the call variants section above.  
vcf_files = "/mnt/grid/gillis/hpc/data/data/sballouz/human/RNAseq/XCGD/pilot/outs/possorted_genomeX_bam.barcoded.cut"


i = 1 
files = read.table(input_files)
vcf  = read.table(vcf_files)


mat.sd = matrix(0, ncol=2, nrow=length(files))
mat.avg = matrix(0, ncol=2, nrow=length(files))

list.x = list()
list.vcf.filts = list() 
list.vcf = list() 
list.skew = list() 
list.skew2 = list() 


i = 1 
for( file in files[,1] ){ 
    #test = read.table(paste0("wigs/",file)) 
    info = file.info(as.character(file))
    if( info$size==0 ) { i = i + 1 ; next;  } 
    test = read.table(file) 
    
    test2 = test[,1:6]
    test2[,2:6] = test2[,2:6] + test[,9:13]
    rm(test)
    colnames(test2) = c( "Pos", "A", "C", "G", "T", "N" )
    sums = rowSums(test2[,2:6])
    b = rowSums(test2[,2:5]>0)
    
    m = match( test2[,1], par1)
    f.p1 = is.na(m)
    m = match( test2[,1], par2)
    f.p2 = is.na(m)
    
    f =  b >= 1 # changed this 
    
    if( sum(f) == 0) { i = i + 1 ; next; }
    
    
    f.all =  f  & f.p1 & f.p2
    d = test2[f.all,2:5]/sums[f.all]
    d[d==0]=NA
    d = as.matrix(d)
    test3.mask = test2[ f.all ,1:5]
   # test3.mask[test3.mask < 2] = 0
    # test3.mask[test3.mask < 10] = 0
    b3.mask = rowSums(test3.mask[,2:5]>0)
    
    if( sum(f.all) == 0) { i = i + 1 ; next; }
    
    m = match(test3.mask[ ,1], vcf[,1])
    f.t = !is.na(m)
    f.v = m[f.t]
    
    if( sum(f.t) < 1) { i = i + 1 ; next; }
    if( sum(f.t) == 1 ){ 
       temp.vcf = cbind(t(as.matrix(d[ , ][f.t]) ) , vcf[f.v,3:4] )
     } else {
      temp.vcf = cbind(d[ ,][f.t,], vcf[f.v,3:4])
     } 
    
    ref1 = as.character(temp.vcf[,5]) == "A"
    ref2 = as.character(temp.vcf[,5]) == "C"
    ref3 = as.character(temp.vcf[,5]) == "G"
    ref4 = as.character(temp.vcf[,5]) == "T"
    
    alt1 = as.character(temp.vcf[,6]) == "A"
    alt2 = as.character(temp.vcf[,6]) == "C"
    alt3 = as.character(temp.vcf[,6]) == "G"
    alt4 = as.character(temp.vcf[,6]) == "T"
    
    ref.only = c(temp.vcf[ref1,1], temp.vcf[ref2,2], temp.vcf[ref3,3], temp.vcf[ref4,4])
    alt.only = c(temp.vcf[alt1,1], temp.vcf[alt2,2], temp.vcf[alt3,3], temp.vcf[alt4,4])
    temp2.vcf = cbind(test3.mask[ ,2:5][f.t,], vcf[f.v,3:4])
    temp3.vcf = temp2.vcf[,1:2]
    temp3.vcf = temp3.vcf* 0
    
    temp3.vcf[ref1,1] = temp2.vcf[ref1,1]
    temp3.vcf[ref2,1] = temp2.vcf[ref2,2]
    temp3.vcf[ref3,1] = temp2.vcf[ref3,3]
    temp3.vcf[ref4,1] = temp2.vcf[ref4,4]
    
    temp3.vcf[alt1,2] = temp2.vcf[alt1,1]
    temp3.vcf[alt2,2] = temp2.vcf[alt2,2]
    temp3.vcf[alt3,2] = temp2.vcf[alt3,3]
    temp3.vcf[alt4,2] = temp2.vcf[alt4,4]
    
    x = c(alt.only, ref.only)
    #	x_f = 2*(1 - ( 0.5 + abs( x - 0.5 )))
    #	y = log(x_f)
    
    # mat.sd[i,1] = sd(x)
    # mat.avg[i,1] = mean(x)
    list.x[[i]] = x
    # mat.sd[i,2] = sd(x)/sqrt(length(x))
    #list.vcf[[i]] = temp3.vcf
    list.vcf[[i]] = cbind(test3.mask[ ,1], temp3.vcf, temp2.vcf)
    i = i + 1 
} 


skip = sapply( 1:length(list.vcf), function(i) length(list.vcf[[i]] ) )  ==0
barcodes = (substr( as.character(files[,1]),6,28 ) )
names(list.vcf) = barcodes
names(list.x) = barcodes

save(list.vcf, list.x, skip, files, barcodes, file="skew_sc.Rdata")




f.x = attr$chr == "chrX"



barcodes = (substr( as.character(files[,1]),6,28 ) )
names(list.vcf) = barcodes
names(list.x) = barcodes

save(list.vcf, list.x, skip, files, barcodes, file="skew_sc.Rdata")

test = do.call(rbind, list.vcf)
col1 = gsub( "\\.[0-9]+", "", rownames(test) )
skewmat = cbind(col1, test)

tests2 = sort(unique(skewmat[,2]))
genes.snps = sapply( 1:length(tests2), function(i) which(tests2[i] >= attr[f.x,2]  &  tests2[i] <= attr[f.x,3] ))

for(i in 1:length(genes.snps)) { if( length(genes.snps[[i]])==0 ) { genes.snps[[i]] = NA } ;  if( length(genes.snps[[i]]) > 1  ) { genes.snps[[i]] = genes.snps[[i]][1]  }   }


genes.snps2 = list() 
for(i in 1:length(genes.snps)) { 
  if( is.na(genes.snps[[i]]) ) { genes.snps2[[i]] = cbind( tests2[i], empty)  } 
  else { 
   genes.snps2[[i]] = cbind( tests2[i],  attr[f.x,][genes.snps[[i]],]  ) 
  } 
}

genes.snps3 = do.call(rbind, genes.snps2)


m2 = grr::matches( skewmat[,2], genes.snps3[,1])  ## this is awesome 
skewmat2 = cbind(skewmat[m2[,1],], genes.snps3[m2[,2],])



save(skewmat2, file="skew_sc_genes.Rdata")



get_haplotype <- function( data){
   pos = data[,1]
   ref = data[,2] > 0
   alt = data[,3] > 0
   hap1 = data[,1] *  0 
   hap2 = hap1
   hap1[ref] = as.character(data[ref,8])
   hap1[alt] = as.character(data[alt,9])
   hap2[ref] = as.character(data[ref,9])
   hap2[alt] = as.character(data[alt,8])
   hap = cbind(pos, hap1, hap2)
   return(hap )
   
}


temp  = skewmat2[,1:2]  
colnames(temp) = c("cell", "bp" ) 
cell_by_snp = make_annotations( temp, unique(temp[,1]), unique(temp[,2] )) 


cell_by_snp.hap = cell_by_snp * NA 

haps = sapply(1:length(skewmat2[,1]), function(i) get_haplotype(skewmat2[i,2:10] ) )

skewmat.merge = cbind( skewmat2,  t(haps[2:3,]))


cell_by_snp.hap = cell_by_snp * NA 
cellids = rownames(cell_by_snp)
snppos = colnames(cell_by_snp) 
colsselect = c(1,2,19,21,22)
for( i in 1:dim(cell_by_snp)[1]){
    subtest = skewmat.merge[skewmat.merge[,1] == cellids[i],colsselect]
    m = match( snppos,subtest[,2] ) 
    f.ss = !is.na(m)
    f.st = m[f.ss]
    cell_by_snp.hap[i,f.ss] = as.character(subtest[f.st,4])
}


temp  = cbind(as.character(skewmat.merge[,1]), paste(skewmat.merge[,2], skewmat.merge[,21] ) ) 
colnames(temp) = c("cell", "bp" ) 
cell_by_hap = make_annotations( temp, unique(temp[,1]), unique(temp[,2] )) 



save(skewmat.merge, cell_by_hap, cell_by_snp.hap, cell_by_snp, file="cell_by.Rdata")
 


bases = c("A","C","G","T")

count_freq <- function(data,list){
 freq = count(data)
 nmax= length(list)
 res = matrix(0, nrow=nmax, ncol=1 )
 rownames(res) = list
 m = match(list, freq[,1])
 f.r = !is.na(m)
 f.f = m[f.r]
 res[f.r,1] = freq[f.f,2]
 return(res)
}



consensus = sapply(1:dim(cell_by_snp.hap)[2],function(i) count_freq(cell_by_snp.hap[,i], bases ) )
rownames(consensus) = bases
colnames(consensus) = colnames(cell_by_snp.hap)

filt = colSums(consensus > 0) == 2 

consensus[,filt]






cell_by_snp.hap.num = cell_by_snp * NA 
cellids = rownames(cell_by_snp)
snppos = colnames(cell_by_snp) 
colsselect = c(1,2,19,21,22)
for( i in 1:dim(cell_by_snp)[1]){
    subtest = skewmat.merge[skewmat.merge[,1] == cellids[i],colsselect]
    m = match( snppos,subtest[,2] ) 
    f.ss = !is.na(m)
    f.st = m[f.ss]
    cell_by_snp.hap.num[i,f.ss] = as.numeric(subtest[f.st,4])
}


cell_by_snp.hap.num[cell_by_snp.hap.num ==1 ] = 0  
cell_by_snp.hap.num[is.na(cell_by_snp.hap.num)] = 0

save(consensus, filt, cell_by_snp.hap.num, file="temp_by_cell.Rdata") 




```
 
