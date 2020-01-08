dir=`pwd`
files=$1
header=header.sam
#for barcode in `cat $files | cut -f2 -d'/' | cut -f1 -d'.' `
#files=../barcodes.counts
#files=../barcodes.tsv

for barcode in `cut -f1 $files `
do
  #echo $barcode
  barcode='CB:Z:'$barcode
  echo $barcode
  cat  $header  $barcode.sam  >  $barcode.2.sam 
  samtools view -bS $barcode.2.sam > $barcode.bam
  samtools index $barcode.bam 
  rm  $barcode.2.sam  
  java -jar /sonas-hs/gillis/hpc/home/sballouz/IGVTools/igvtools.jar  count -z 0 -w 1 --bases --strands read $barcode.bam  $barcode.wig /sonas-hs/gillis/hpc/home/sballouz/lyon/XSKEW/sc/10x/v3_chemistry/10k_pbmcs_healthy_donor/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa 
  mv  $barcode.wig ../wigs/
done
