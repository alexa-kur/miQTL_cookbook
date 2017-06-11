#!/bin/bash
trap 'kill $(jobs -p)' SIGINT 
echo "First step.."
date

mkdir tmp

InFile=$1
cutoff=$2

cat $InFile |sed 's/-//' |tr -s '\t' '\t'|sed 's/ //g' |sed -e 's/\[/./g;s/\]//g'|perl -pe 's/[.][.]/./g' > tmp/tmp1.out

echo "Second step.."
date
awk -F"\t" '{print $1"\t"$3"."$2"\t"$4"\t"$6"."$5"\t"$7"\t"$9"."$8"\t"$10"\t"$12"."$11"\t"$13"\t"$15"."$14"\t"$16"\t"$18"."$17"\t"$19"\t"$21"."$20"\t"$22}' tmp/tmp1.out > tmp/tmp2.out

echo "Third step.."
date
awk -F $'\t' -v cutoff="$cutoff" 'BEGIN {OFS = FS}{
	if($3 < cutoff) {$2 = "NOTAX_root"}
	if($5 < cutoff) {$4 = "NOTAX_domain"}
	if($7 < cutoff) {$6 = "NOTAX_phylum"}
	if($9 < cutoff) {$8 = "NOTAX_class"}
	if($11 < cutoff) {$10 = "NOTAX_order"}
	if($13 < cutoff) {$12 = "NOTAX_family"}
	if($15 < cutoff) {$14 = "NOTAX_genus"}
	print $0 }' tmp/tmp2.out |cut -f-2,4,6,8,10,12,14 > tmp/tmp3.out
echo "Fourth step.."
date

for i in {2..8}; do for j in $(cat tmp/tmp3.out |awk -F"\t" -v i="$i" '{print $i}' |sort |uniq );do printf SampleID"\t"$j"\n" > tmp/$j.'txt'|grep -w $j tmp/tmp3.out |awk '{print $1}' |sed 's/_/ /' |awk '{print $1}' |sort |uniq -c |awk '{print $2"\t"$1}' >> tmp/$j.'txt'; done; done
date

cd tmp
echo -e "filenames=list.files(pattern='.txt', full.names=TRUE)\ndatalist = lapply(filenames, function(x){read.table(file=x,header=T, sep='\t')})\ndata = Reduce(function(x,y) {merge(x,y, by='SampleID', all=TRUE)}, datalist)\ndata[is.na(data)] <- 0\nwrite.table(data, 'taxonomy_table.txt', col.names=T, row.names=F, quote=F, sep='\t')" > myscript.R

echo "Fifth step.. <Making table>"
Rscript myscript.R

date

mv taxonomy_table.txt ../.

echo "file taxonomy_table.txt has been made!"

cd ..
rm -rf tmp
