#! /usr/bin/perl


#@BARCODES = "cat /sonas-hs/gillis/hpc/home/sballouz/lyon/XSKEW/sc/10x/v3_chemistry/10k_pbmcs_healthy_donor/filtered_feature_bc_matrix/barcodes.tsv";

#@HASH{@BARCODES} = 1;
# possorted_genome_bam_sorted_chrX.split.filtered.bam 
$IN=$ARGV[0];


open (IN, "samtools view $IN | ");
while(<IN>){
  chomp $_ ;
  @temp = split(/\t/,$_);
  $barcodeID = $temp[11];
  print "$barcodeID\n";
  push(@{$BARCODES{$barcodeID}}, $_);

}
close IN; 

#$header = 'cat header.sam'; 

my @ke = keys %BARCODES;

foreach (@ke) {
  $out = $_;
  open (OUT, "> $out.sam");
  @temp = @{$BARCODES{$_}};
#  print OUT "$header";
  print OUT join "\n", @temp;
  close OUT;
}
