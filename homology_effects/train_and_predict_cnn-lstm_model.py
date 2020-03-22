#!/usr/bin/env python2.7
'''
amp_train_and_predict_cnn-lstm_model.py By: Dan Veltri
Takes a protein multi-FASTA file as input and outputs prediction values for each.
Assumes peptides are >= 10 and <= 200 AA in length (AA's longer than 'max_length' are ignored).

Prediction probabilities >0.5 signifiy AMPs, <= 0.5 signify non-AMPs.

User should provide training, validation, and testing FASTA files for AMPs and decoys, respectivley. 
Ensure 'max_length' is >= to the longest peptide in those files.

NOTE on a multi-threaded machine TF will still produce stochastic results even if you set a random seed.
'''
from __future__ import print_function # enable Python3 printing
import numpy as np
from sklearn.metrics import confusion_matrix, roc_auc_score, matthews_corrcoef, classification_report
from sklearn.utils import shuffle
from Bio import SeqIO
import tensorflow as tf

from keras.models import Sequential, load_model, model_from_json
from keras.layers import Dense, LSTM
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.layers.embeddings import Embedding
from keras.preprocessing import sequence


#User-set variables
max_length = 200
amp_train_fasta       = 'AMP.tr.fa'
amp_validate_fasta    = 'AMP.eval.fa'
amp_test_fasta        = 'AMP.te.fa'
decoy_train_fasta     = 'DECOY.tr.fa'
decoy_validate_fasta  = 'DECOY.eval.fa'
decoy_test_fasta      = 'decoy.te.fa'

#Model params
embedding_vector_length = 128
nbf = 64 		# No. Conv Filters
flen = 16 		# Conv Filter length 
nlstm = 100 	# No. LSTM layers
ndrop = 0.1     # LSTM layer dropout
nbatch = 32 	# Fit batch No.
nepochs = 10    # No. training rounds

amino_acids = "XACDEFGHIKLMNPQRSTVWY"
aa2int = dict((c, i) for i, c in enumerate(amino_acids))

X_train = []
y_train = []
X_val = []
y_val = []
X_test = []
y_test = []

print("Encoding training/testing sequences...")
for s in SeqIO.parse(amp_train_fasta,"fasta"):
    X_train.append([aa2int[aa] for aa in str(s.seq).upper()])
    y_train.append(1)
for s in SeqIO.parse(amp_validate_fasta,"fasta"):
    X_val.append([aa2int[aa] for aa in str(s.seq).upper()])
    y_val.append(1)
for s in SeqIO.parse(amp_test_fasta,"fasta"):
    X_test.append([aa2int[aa] for aa in str(s.seq).upper()])
    y_test.append(1)
for s in SeqIO.parse(decoy_train_fasta,"fasta"):
    X_train.append([aa2int[aa] for aa in str(s.seq).upper()])
    y_train.append(0)
for s in SeqIO.parse(decoy_validate_fasta,"fasta"):
    X_val.append([aa2int[aa] for aa in str(s.seq).upper()])
    y_val.append(0)  
for s in SeqIO.parse(decoy_test_fasta,"fasta"):
    X_test.append([aa2int[aa] for aa in str(s.seq).upper()])
    y_test.append(0)

# Pad input sequences
X_train = sequence.pad_sequences(X_train, maxlen=max_length)
X_val = sequence.pad_sequences(X_val, maxlen=max_length)
X_test = sequence.pad_sequences(X_test, maxlen=max_length)

# Shuffle training sequences
X_train, y_train = shuffle(X_train, np.array(y_train))
X_val, y_val = shuffle(X_val, np.array(y_val))

print("Compiling model...")
model = Sequential()
model.add(Embedding(21, embedding_vector_length, input_length=max_length))
model.add(Conv1D(filters=nbf, kernel_size=flen, padding="same", activation='relu'))
model.add(MaxPooling1D(pool_size=5))
model.add(LSTM(nlstm, use_bias=True, dropout=ndrop, return_sequences=False))#,merge_mode='ave'))
model.add(Dense(1, activation='sigmoid'))
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

print("Training now...")
model.fit(X_train, np.array(y_train), epochs=nepochs, batch_size=nbatch, verbose=1)

print("\nGathering Testing Results...")
preds = model.predict(X_test)
pred_class = np.rint(preds) #round up or down at 0.5
true_class = np.array(y_test)
tn, fp, fn, tp = confusion_matrix(true_class,pred_class).ravel()
roc = roc_auc_score(true_class,preds) * 100.0
mcc = matthews_corrcoef(true_class,pred_class)
acc = (tp + tn) / (tn + fp + fn + tp + 0.0) * 100.0
sens = tp / (tp + fn + 0.0) * 100.0
spec = tn / (tn + fp + 0.0) * 100.0
prec = tp / (tp + fp + 0.0) * 100.0

print("\nTP\tTN\tFP\tFN\tSens\tSpec\tAcc\tMCC\tauROC\tPrec")
print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(tp,tn,fp,fn,np.round(sens,4),np.round(spec,4),np.round(acc,4),np.round(mcc,4),np.round(roc,4),np.round(prec,4)))

# END PROGRAM
