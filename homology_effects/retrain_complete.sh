# Clustering databases
zcat old_files/*.faa.gz > tmp
cdhit -i tmp -o tmp_clusters -c 0.8 -b 10 -T 3 -l 9 -d 0 -aS 0.9 -p 1 -g 0 -t0
rm -rf tmp

###### Testing clustering efficiency
################################################################################################################################################################
# Formatting blastdb
makeblastdb -in tmp_clusters -dbtype prot -out reference

# Smith-waterman alignment search
swipe -d reference -i tmp_clusters -a 3 -m '8 std qcovs' -o sw_swipe -p 1
mkdir clustering
mv tmp_clusters* clustering/

# Listing names
grep ">" clustering/tmp_clusters | awk '{print "gnl|BL_ORD_ID|"NR-1"\t"$1}' | sed 's/>//g' > ref.names.list
grep -v ">" clustering/tmp_clusters | awk '{print length}' > t
paste -d'\t' ref.names.list t > t2; rm -rf t; mv t2 ref.names.list

# Parsing output
awk '$3 >= 80 && $11 <= 1e-5' sw_swipe > tmp
# Fixing names
while read a b; do sed -i "s/$a\t/$b\t/g" tmp; done < ref.names.list
awk '$1 != $2' tmp | sort -k1,1 > tmp.2
sort -k2,2 ref.names.list | cut -f2,3 > ref
join ref tmp.2 | sed 's/ /\t/g' | awk '{print $0"\t"$6/$2"\t"$6/$4}' > tmp.3

# Filtering of seqs with less than 90% of the shorter sequence
awk '$15 >= 0.9' tmp.3 > tmp.3.filt
awk '$16 >= 0.9' tmp.3 > tmp.3.filt2

# Finishing process
cat tmp.3.* | sort -k1,1 | uniq > sw_swipe.parsed; rm -rf tmp*
## final parsed table has the following fields: 
#query q_len subject s_length %id alignment_length mismatches gap_openings query_start query_end subject_start subject_end E_value bit_score s_covs q_covs
mkdir SW
mv sw_* SW
rm -rf *pin *phr *psq tmp reference ref*
################################################################################################################################################################

############ Testing models in MACREL
################################################################################################################################################################
# Creating databases
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < clustering/tmp_clusters | sort -k1,1 > tmp.1

grep ">NAMP_" tmp.1 > namp
grep -v ">NAMP_" tmp.1 > amp
rm -rf tmp.1
wc -l *amp

#######################################################
##### Check results                               #####
##  1697 AMPs                                     #####
## 90956 NAMPs                                    #####
## Then we divided data sets like:                #####
##   500 Test seqs (1:1)                          #####
## Train sets with all AMPs except those in train #####
## Train sets with different ammounts of NAMPs    #####
#######################################################

shuf -n 500 namp > namptest
shuf -n 500 amp > amptest
cat amptest namptest > test_set_1_1
grep -v -f namptest namp > namptrain
grep -v -f amptest amp > amptrain
shuf -n 1197 namptrain > t
cat t amptrain > 1_1_trainset
shuf -n 6000 namptrain > t
cat t amptrain > 1_5_trainset
shuf -n 12000 namptrain > t
cat t amptrain > 1_10_trainset
shuf -n 24000 namptrain > t
cat t amptrain > 1_20_trainset
shuf -n 36000 namptrain > t
cat t amptrain > 1_30_trainset
shuf -n 48000 namptrain > t
cat t amptrain > 1_40_trainset
shuf -n 60000 namptrain > t
cat t amptrain > 1_50_trainset
rm -rf t amptrain namptrain namptest amptest

# calculating features
for i in *_trainset; do awk '{print $1"\n"$2}' $i > t; mv t $i; python3 AMP_features.py $i $i.tsv; done
awk '{print $1"\n"$2}' test_set_1_1 > t; mv t test_set_1_1; python3 AMP_features.py test_set_1_1 test_set_1_1.tsv
mkdir sets
mv *_trainset sets
mv test_set_1_1 sets

# You should fix manually the tsv files in the field "group"
# Then, keep the following:
sed -i 's/,/\./g' *.tsv

# Training models and calculating confusion matrices
python3 train-models.py
mkdir features
mv *tsv features

# homology classification
mkdir homolog_class
cd homolog_class
for i in ../sets/*_trainset; do makeblastdb -in $i -dbtype prot -out ${i}; done
for i in ../sets/*_trainset; do  blastp -db ../sets/${i} -query ../sets/test_set_1_1 -out homolog_class -evalue 1e-5 -word_size 5 -qcov_hsp_perc 90 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -window_size 10 -subject_besthit -max_target_seqs 1 -num_threads 3; name=`echo $i | sed 's/\.\.\/sets\///g'`; mv homolog_class homolog_class.$name; done
for i in homolog_class*; do awk '$3 >= 50 && $11 <= 1e-5 && $13 >= 80' $i > $i.parsed; done
for i in *.parsed; do cut -f1,2 $i | sed 's/_/\t/g' | awk '{print $1"_"$3}' | sort | uniq -c > ${i/.parsed/.calc}; done
# to determine parameters of performance in classification
# TP were given as "AMP_AMP"
# TN were given as "NAMP_NAMP"
# FP were given as (500 -  TN)
# FN were given as (500 - TP)
for i in *.calc; do k=500; TP=`grep -w "AMP_AMP" $i | awk '{print $1}'`; FN=`expr $k - $TP`; TN=`grep -w "NAMP_NAMP" $i | awk '{print $1}'`; FP=`expr $k - $TN`; echo -e "${i/.calc/}\t$TP\t$FP\t$TN\t$FN"; done
rm -rf *parsed *calc
cd ../

### Training iAMP2L
## Calculate features
R --vanilla --slave < iAMP2L.features.R
## Convert tables to feature tables
grep ">" sets/test_set_1_1 | sed 's/_.*//g' | sed 's/>//g' | sed 's/ //g' > t; paste -d'\t' t test_set_1_1.tmp > test_set_1_1.tsv; rm -rf test_set_1_1.tmp t
for i in $(ls sets/*_trainset); do name=`echo $i | sed 's/sets\///g'`; grep ">" $i | sed 's/_.*//g' | sed 's/>//g' | sed 's/ //g' > t; paste -d'\t' t $name.tmp > $name.tsv; rm -rf $i.tmp t; done
## training models
# dividing datasets
mkdir feat_iAMP2L
mkdir sepdata
for i in $(ls *.tsv); do cut -f1 $i > sepdata/y.$i; cut -f2- $i > sepdata/x.$i; mv $i feat_iAMP2L/; done
# training
python train_iAMP2L.py
# cleaning
rm -rf *.tmp

### Training AMP Scanner
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < sets/test_set_1_1 | grep ">AMP_" > AMP.te.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < sets/test_set_1_1 | grep ">NAMP_" > decoy.te.fa

grep -v -w -f AMP.te.fa amp > amp.tr
shuf -n 600 amp.tr > AMP.tr.fa
grep -v -w -f AMP.tr.fa amp.tr > AMP.eval.fa

grep -v -w -f decoy.te.fa namp > namp.tr
shuf -n 45228 namp.tr > decoy.tr.fa
grep -v -w -f decoy.tr.fa namp.tr > decoy.eval.fa

for i in $(ls *.fa); do awk '{print $1"\n"$2}' $i > tm; mv tm $i; done
rm -rf amp.tr namp.tr amp namp

### Preparing data for testing 1:1
head -600 decoy.eval.fa > DECOY.eval.fa 
head -600 decoy.tr.fa > DECOY.tr.fa
python train_and_predict_cnn-lstm_model.py > test_1_1

### Preparing data for testing 1:5
head -3000 decoy.eval.fa > DECOY.eval.fa 
head -3000 decoy.tr.fa > DECOY.tr.fa
python train_and_predict_cnn-lstm_model.py > test_1_5

### Preparing data for testing 1:10
head -6000 decoy.eval.fa > DECOY.eval.fa 
head -6000 decoy.tr.fa > DECOY.tr.fa
python train_and_predict_cnn-lstm_model.py > test_1_10

### Preparing data for testing 1:20
head -12000 decoy.eval.fa > DECOY.eval.fa 
head -12000 decoy.tr.fa > DECOY.tr.fa
python train_and_predict_cnn-lstm_model.py > test_1_20

### Preparing data for testing 1:30
head -18000 decoy.eval.fa > DECOY.eval.fa 
head -18000 decoy.tr.fa > DECOY.tr.fa
python train_and_predict_cnn-lstm_model.py > test_1_30

### Preparing data for testing 1:40
head -24000 decoy.eval.fa > DECOY.eval.fa 
head -24000 decoy.tr.fa > DECOY.tr.fa
python train_and_predict_cnn-lstm_model.py > test_1_40

### Preparing data for testing 1:50
head -30000 decoy.eval.fa > DECOY.eval.fa 
head -30000 decoy.tr.fa > DECOY.tr.fa
python train_and_predict_cnn-lstm_model.py > test_1_50

mkdir AMP_Scanner
mv AMP.* AMP_Scanner/
mv decoy.* AMP_Scanner/
mv test_1_* AMP_Scanner/

rm -rf DECOY.* model.*
