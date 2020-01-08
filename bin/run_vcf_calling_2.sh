picard=/sonas-hs/gillis/hpc/home/sballouz/GenomeAnalysisTK/picard.jar
GATK=/sonas-hs/gillis/hpc/home/sballouz/GenomeAnalysisTK/GenomeAnalysisTK.jar
igvtools=/sonas-hs/gillis/hpc/home/sballouz/IGVTools/igvtools.jar
samtools=/opt/hpc/bin/samtools
genome_fa=/sonas-hs/gillis/hpc/home/sballouz/STAR/GRCh38_Gencode25/GRCh38.p7.genome.fa


inDir=$1
chr=$2
outPrefix=`pwd`
echo $outPrefix
inI=$((SGE_TASK_ID-1))
echo $inI
outDir=$inDir
cd $outDir

echo "Adding read groups to bam file"
java -jar /sonas-hs/gillis/hpc/home/sballouz/GenomeAnalysisTK/picard.jar AddOrReplaceReadGroups \
    I=possorted_genome_bamX.filt.bam \
    O=possorted_genome_bamX.rg.bam \
    SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

echo "Marking duplicates"
java -jar /sonas-hs/gillis/hpc/home/sballouz/GenomeAnalysisTK/picard.jar MarkDuplicates \
    I=possorted_genome_bamX.rg.bam  \
    O=possorted_genome_bamX.dedupped.bam \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

echo "Splitting and trimming"
java -jar /sonas-hs/gillis/hpc/home/sballouz/GenomeAnalysisTK/GenomeAnalysisTK.jar \
   -T SplitNCigarReads         \
   -R /sonas-hs/gillis/hpc/home/sballouz/lyon/XSKEW/sc/10x/v3_chemistry/10k_pbmcs_healthy_donor/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
   -I possorted_genome_bamX.dedupped.bam \
   -o possorted_genome_bamX.split.filtered.bam \
   -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS


echo "Haplotype calling (?)"
java -jar /sonas-hs/gillis/hpc/home/sballouz/GenomeAnalysisTK/GenomeAnalysisTK.jar \
   -T HaplotypeCaller \
   -R /sonas-hs/gillis/hpc/home/sballouz/lyon/XSKEW/sc/10x/v3_chemistry/10k_pbmcs_healthy_donor/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
   -L X \
   -I possorted_genome_bamX.split.filtered.bam \
   -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 \
   -o possorted_genome_bamX.split.filtered.vcf

echo "Counting SNPs"
java -jar /sonas-hs/gillis/hpc/home/sballouz/IGVTools/igvtools.jar count -z 0 -w 1 --bases --strands read   \
   possorted_genome_bamX.split.filtered.bam   \
   possorted_genome_bamX.split.filtered.wig   \
   /sonas-hs/gillis/hpc/home/sballouz/lyon/XSKEW/sc/10x/v3_chemistry/10k_pbmcs_healthy_donor/refdata-cellranger-GRCh38-3.0.0/fasta/X.fa

 

	grep ^X possorted_genome_bamX.split.filtered.vcf  | cut -f2-5 > A
	cut -f1 A > B
	grep -Fw -f B /sonas-hs/gillis/hpc/home/sballouz/lyon/ncbi/dbGAP-11987/C > D
	grep -Fw -f D possorted_genome_bamX.split.filtered.wig  > possorted_genome_bamX.intersect.wig 
	grep -Fw -f D A > possorted_genome_bamX.cut
	rm A B D 




