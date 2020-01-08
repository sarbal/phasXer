dir=`pwd`
chr=chrX
files=$1
for file in `cat $files | cut -f1`
do
	echo $file $chr	
	cd $dir/$file
	bash ../bin/run_vcf_calling_2.sh $dir/$file/ 
	cd ../
done
