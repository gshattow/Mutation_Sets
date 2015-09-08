#!/bin/bash   

# This script runs on the output of the VEP which includes the sift scores.
# For each chromosome, it moves and unzips the files.
# Then taking only the lines with sift scores, it records the line number 
# (minus the header) it came from.
# It looks up the ENSEMBL and HGNC gene cataplgues corresponding to each rs number
# using symbol_lookup.py. (optional)
# Then it combines all the columns into 
# line number, rs number, ENSEMBL ID, (HGNC ID), sift score, and phenotype score.

dataset='20130502'
dataloc='ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/'$dataset'/supporting/functional_annotation/filtered/'
dir='../data/F'$dataset'/'
mkdir $dir
time for c in {1..22}
do
	file='ALL.chr'$c'.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.annotation.vcf.gz'
	wget -P $dir $dataloc$file
	newfile=$dir'ALL.chr'$c'.vcf.gz'
	echo 'moving '$dir$file'to '$newfile
	mv $dir$file $newfile
	echo 'unzipping '$newfile
	gunzip $newfile
	unzipped=${newfile::-3}
	siftfile=$dir'SIFT.chr'$c'.vcf'
	header=$(head -n 1000 $unzipped | grep "#" | wc -l)
	echo 'taking columns from '$unzipped
	grep -n "deleterious\|tolerated" $unzipped | grep "rs" > $siftfile
	lines=$dir'line.txt'
	ids=$dir'id.txt'
	info=$dir'info.txt'
	sifts=$dir'sift.txt'
	awk '{print $1}' $siftfile | awk -F ":" '{print $1-'$header'}' > $lines #.txt
	awk '{print $3}' $siftfile > $ids #.txt
	awk '{print $8}' $siftfile > $info #.txt
	awk -F "|" '{print $5"\t"$17"\t"$18}' $info |\
		sed 's/(/\t/g' | sed 's/)//g' > $sifts

###### Commented out lines are for symbol lookups to shift from ENSG
###### numbers to HGNC, etc. Which is now unnecessary, but I may want
###### to work it back in at some point
#	cat ../symbol_lookup.py | sed "s/XXX/$c/g" > sym_lookup.py
#	python sym_lookup.py

	final=$dir'sifted.SIFT.chr'$c'.txt'
#	pr -m -t -s\  line.txt id.txt sift.txt hgnc_symbols.txt | gawk '{print $1,$2,$3,$8,$5,$7}' > $final
	pr -m -t -s\  $lines $ids $sifts | gawk '{print $1,$2,$3,$5,$7}' > $final
	
#	echo 'line, id, ENSG id, HGNC sym, SIFT, and phen printed to '$final
	echo 'line, id, ENSG id, SIFT, and phenotype printed to '$final
	rm $unzipped
	rm $siftfile
	rm $lines
	rm $ids
	rm $info
	rm $sifts
	
done

plotsdir='../plots'
mkdir $plotsdir