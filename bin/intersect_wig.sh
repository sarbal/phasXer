dir=`pwd`
files=$1 
#barcodes.tsv
VCF=$2
#../possorted_genome_bam_sorted_chrX.filtered.vcf
SNPs=/sonas-hs/gillis/hpc/home/sballouz/lyon/ncbi/dbGAP-11987/C
cut -f2-5 $VCF > A
cut -f1 A > B
grep -Fw -f B $SNPs  > D
#grep -Fw -f D A > possorted_genomeX_bam.barcoded.cut

for barcode in `cut -f1 $files `
do
  barcode='CB:Z:'$barcode
  echo $barcode
  grep -Fw -f D wigs/$barcode.wig > wigs/$barcode.intersect.wig 
done
