#!/usr/bin/env bash

set -e

source ~/anaconda3/etc/profile.d/conda.sh
conda activate

BENCHMARK_DIR=$PWD/BENCHMARK
spurio="/path/to/spurio.py"
db="/path/to/spurio/db/fullsource_filter.fa"

##################################################################################################################################################################
# Creating directories

mkdir -p $BENCHMARK_DIR/metagenomes/{data macrel}
mkdir -p $BENCHMARK_DIR/metagenomes/data/{metaT bam profiles}
mkdir -p $BENCHMARK_DIR/metagenomes/macrel/{BASHLOG IDS AMP LOG CONTIGS}
out="$BENCHMARK_DIR/metagenomes/macrel/"

##################################################################################################################################################################
# To download metagenomes and metatranscriptomes

for i in $(cat SRR_Acc_List.txt); 
do 
	fastq-dump -I --split-files --outdir $BENCHMARK_DIR/metagenomes/data/ --gzip $i
done

for i in $(cat SRR_Acc_List2.txt); 
do 
	fastq-dump -I --split-files --outdir $BENCHMARK_DIR/metagenomes/data/metaT/ --gzip $i
done

##################################################################################################################################################################
# Processing files with macrel
while read a
do

	echo -e "Doing $a metagenome"
	$macrel reads -1 $BENCHMARK_DIR/metagenomes/data/"$a"_1.fastq.gz -2 $BENCHMARK_DIR/metagenomes/data/"$a"_2.fastq.gz 
	\--outfolder $out --outtag $a > $out/bashlog.$a.txt
	
	mv $a.log $out/LOG/
	mv $a.ids.tsv.gz $out/IDS/
	mv $a.tsv.gz $out/AMP/
	mv $a.fna.gz $out/CONTIGS/
	
done < SRR_Acc_List.txt

##################################################################################################################################################################
# Then all AMPs were pooled:

cd $out/AMP/

zcat * | sed '/Access/d' | cut -f2- | sort -k1,1 | uniq > tmp

echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability" > header

cat header tmp | pigz --best > heinz2016.AMP.tsv.gz

rm -rf header tmp

##################################################################################################################################################################
# Spurious analysis - To this it was used the tool from Hops et al. (2018), more info in: <https://bitbucket.org/bateman-group/spurio/src/master/>

zcat heinz2016.AMP.tsv.gz | sed '1,1d' | awk '{print ">"NR"\n"$1}' > tmp.fa
python3 $spurio -s 1 -e $(grep -c ">" tmp.fa) -v 1 -r $db -q tmp.fa -qt spurio_res
awk '$2 > 0.8' spurious_res.txt | awk '{print $1}' > list
zgrep -v -w -A1 -f list tmp.fa | grep -v ">" > AMP
rm -rf list tmp.fa
zgrep -w -f AMP heinz2016.AMP.tsv.gz > heinz2016.AMP.tsv
rm -rf AMP heinz2016.AMP.tsv.gz
pigz --best heinz2016.AMP.tsv

##################################################################################################################################################################
##################################################################################################################################################################
# Abundance

while read a
do
	echo -e "Doing $a metagenome"
	$macrel abundances --fasta tmp.fa -1 $BENCHMARK_DIR/metagenomes/data/metaT/$a_1.fastq.gz -2 $BENCHMARK_DIR/metagenomes/data/metaT/$a_2.fastq.gz 
	\--outfolder $BENCHMARK_DIR/metagenomes/data/bam/ --outtag $a > $BENCHMARK_DIR/metagenomes/data/bam/bashlog.$a.txt

done < SRR_Acc_List2.txt
