from __future__ import print_function # enable Python3 printing

import pprint
import operator

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf

from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.datasets import load_iris, load_breast_cancer
from sklearn.metrics import confusion_matrix, roc_auc_score, matthews_corrcoef, classification_report, accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils import shuffle
from keras.models import Sequential, load_model, model_from_json

#User-set variables
xtest            = pd.read_csv('sepdata/x.test_set_1_1.tsv', sep='\t', index_col=0).values
x11train         = pd.read_csv('sepdata/x.1_1_trainset.tsv', sep='\t', index_col=0).values
x15train         = pd.read_csv('sepdata/x.1_5_trainset.tsv', sep='\t', index_col=0).values
x110train        = pd.read_csv('sepdata/x.1_10_trainset.tsv', sep='\t', index_col=0).values
x120train        = pd.read_csv('sepdata/x.1_20_trainset.tsv', sep='\t', index_col=0).values
x130train        = pd.read_csv('sepdata/x.1_30_trainset.tsv', sep='\t', index_col=0).values
x140train        = pd.read_csv('sepdata/x.1_40_trainset.tsv', sep='\t', index_col=0).values
x150train        = pd.read_csv('sepdata/x.1_50_trainset.tsv', sep='\t', index_col=0).values
ytest            = pd.read_csv('sepdata/y.test_set_1_1.tsv', sep='\t', index_col=0)
ytest1           = ytest.index.values.tolist()
y11train         = pd.read_csv('sepdata/y.1_1_trainset.tsv', sep='\t', index_col=0)
y11train1        = y11train.index.values.tolist()
y15train         = pd.read_csv('sepdata/y.1_5_trainset.tsv', sep='\t', index_col=0)
y15train1        = y15train.index.values.tolist()
y110train        = pd.read_csv('sepdata/y.1_10_trainset.tsv', sep='\t', index_col=0)
y110train1       = y110train.index.values.tolist()
y120train        = pd.read_csv('sepdata/y.1_20_trainset.tsv', sep='\t', index_col=0)
y120train1        = y120train.index.values.tolist()
y130train        = pd.read_csv('sepdata/y.1_30_trainset.tsv', sep='\t', index_col=0)
y130train1        = y130train.index.values.tolist()
y140train        = pd.read_csv('sepdata/y.1_40_trainset.tsv', sep='\t', index_col=0)
y140train1        = y140train.index.values.tolist()
y150train        = pd.read_csv('sepdata/y.1_50_trainset.tsv', sep='\t', index_col=0)
y150train1        = y150train.index.values.tolist()

#Defining FKNN (https://github.com/sahilsehwag/FuzzyKNN/blob/master/fuzzy_knn.ipynb)
class FuzzyKNN(BaseEstimator, ClassifierMixin):
    def __init__(self, k=3, plot=False):
        self.k = k
        self.plot = plot
        
        
    def fit(self, X, y=None):
        self._check_params(X,y)
        self.X = X
        self.y = y
        
        self.xdim = len(self.X[0])
        self.n = len(y)
        
        classes = list(set(y))
        classes.sort()
        self.classes = classes
        
        self.df = pd.DataFrame(self.X)
        self.df['y'] = self.y
        
        self.memberships = self._compute_memberships()
        
        self.df['membership'] = self.memberships
        
        self.fitted_ = True
        return self
    
    
    def predict(self, X):
        if self.fitted_ == None:
            raise Exception('predict() called before fit()')
        else:
            m = 2
            y_pred = []
            
            for x in X:
                neighbors = self._find_k_nearest_neighbors(pd.DataFrame.copy(self.df), x)
                
                votes = {}
                for c in self.classes:
                    den = 0
                    for n in range(self.k):
                        dist = np.linalg.norm(x - neighbors.iloc[n,0:self.xdim])
                        den += 1 / (dist ** (2 / (m-1)))
                    
                    neighbors_votes = []
                    for n in range(self.k):
                        dist = np.linalg.norm(x - neighbors.iloc[n,0:self.xdim])
                        num = (neighbors.iloc[n].membership[c]) / (dist ** (2 / (m-1)))
                        
                        vote = num/den
                        neighbors_votes.append(vote)
                    votes[c] = np.sum(neighbors_votes)
                    
                pred = max(votes.items(), key=operator.itemgetter(1))[0]
                y_pred.append((pred, votes))
                
            return y_pred
        
        
    def score(self, X, y):
        if self.fitted_ == None:
            raise Exception('score() called before fit()')
        else:
            predictions = self.predict(X)
            y_pred = [t[0] for t in predictions]
            confidences = [t[1] for t in predictions]
            
            return accuracy_score(y_pred=y_pred, y_true=y)
    
        
    def _find_k_nearest_neighbors(self, df, x):
        X = df.iloc[:,0:self.xdim].values
        
        df['distances'] = [np.linalg.norm(X[i] - x) for i in range(self.n)]
        
        df.sort_values(by='distances', ascending=True, inplace=True)
        neighbors = df.iloc[0:self.k]
        
        return neighbors

                
    def _get_counts(self, neighbors):
        groups = neighbors.groupby('y')
        counts = {group[1]['y'].iloc[0]:group[1].count()[0] for group in groups}
        
        return counts
        
        
    def _compute_memberships(self):
        memberships = []
        for i in range(self.n):
            x = self.X[i]
            y = self.y[i]
            
            neighbors = self._find_k_nearest_neighbors(pd.DataFrame.copy(self.df), x)
            counts = self._get_counts(neighbors)
        
            membership = dict()
            for c in self.classes:
                try:
                    uci = 0.49 * (counts[c] / self.k)
                    if c == y:
                        uci += 0.51
                    membership[c] = uci
                except:
                    membership[c] = 0
                    
            memberships.append(membership)
        return memberships
        
        
    def _check_params(self, X, y):
        if type(self.k) != int:
            raise Exception('"k" should have type int')
        elif self.k >= len(y):
            raise Exception('"k" should be less than no of feature sets')
        elif self.k % 2 == 0:
            raise Exception('"k" should be odd')
            
        if type(self.plot) != bool:
            raise Exception('"plot" should have type bool')

custModel = FuzzyKNN()
ytest2 = [1 if i == 'AMP' else 0 for i in ytest1]
true_class = np.array(ytest2)

#### Proportion 1:1
X = x11train
y = [1 if i == 'AMP' else 0 for i in y11train1]
xTrain, xTest, yTrain, yTest = train_test_split(X,y)
custModel.fit(xTrain, yTrain)

print("\nGathering Testing Results...")
preds = custModel.predict(xtest)
classe = [pre_class[0] for pre_class in preds]
pred_class = np.array(classe)
tn, fp, fn, tp = confusion_matrix(true_class,pred_class).ravel()

print("\nTP\tFP\tTN\tFN")
print("{}\t{}\t{}\t{}".format(tp,fp,tn,fn))

#### Proportion 1:5

X = x15train
y = [1 if i == 'AMP' else 0 for i in y15train1]
xTrain, xTest, yTrain, yTest = train_test_split(X,y)
custModel.fit(xTrain, yTrain)

preds = custModel.predict(xtest)
classe = [pre_class[0] for pre_class in preds]
pred_class = np.array(classe)
tn, fp, fn, tp = confusion_matrix(true_class,pred_class).ravel()
print("{}\t{}\t{}\t{}".format(tp,fp,tn,fn))

#### Proportion 1:10

X = x110train
y = [1 if i == 'AMP' else 0 for i in y110train1]
xTrain, xTest, yTrain, yTest = train_test_split(X,y)
custModel.fit(xTrain, yTrain)

preds = custModel.predict(xtest)
classe = [pre_class[0] for pre_class in preds]
pred_class = np.array(classe)
tn, fp, fn, tp = confusion_matrix(true_class,pred_class).ravel()
print("{}\t{}\t{}\t{}".format(tp,fp,tn,fn))

#### Proportion 1:20

X = x120train
y = [1 if i == 'AMP' else 0 for i in y120train1]
xTrain, xTest, yTrain, yTest = train_test_split(X,y)
custModel.fit(xTrain, yTrain)

preds = custModel.predict(xtest)
classe = [pre_class[0] for pre_class in preds]
pred_class = np.array(classe)
tn, fp, fn, tp = confusion_matrix(true_class,pred_class).ravel()
print("{}\t{}\t{}\t{}".format(tp,fp,tn,fn))

#### Proportion 1:30

X = x130train
y = [1 if i == 'AMP' else 0 for i in y130train1]
xTrain, xTest, yTrain, yTest = train_test_split(X,y)
custModel.fit(xTrain, yTrain)

preds = custModel.predict(xtest)
classe = [pre_class[0] for pre_class in preds]
pred_class = np.array(classe)
tn, fp, fn, tp = confusion_matrix(true_class,pred_class).ravel()
print("{}\t{}\t{}\t{}".format(tp,fp,tn,fn))

#### Proportion 1:40

X = x140train
y = [1 if i == 'AMP' else 0 for i in y140train1]
xTrain, xTest, yTrain, yTest = train_test_split(X,y)
custModel.fit(xTrain, yTrain)

preds = custModel.predict(xtest)
classe = [pre_class[0] for pre_class in preds]
pred_class = np.array(classe)
tn, fp, fn, tp = confusion_matrix(true_class,pred_class).ravel()
print("{}\t{}\t{}\t{}".format(tp,fp,tn,fn))

#### Proportion 1:50

X = x150train
y = [1 if i == 'AMP' else 0 for i in y150train1]
xTrain, xTest, yTrain, yTest = train_test_split(X,y)
custModel.fit(xTrain, yTrain)

preds = custModel.predict(xtest)
classe = [pre_class[0] for pre_class in preds]
pred_class = np.array(classe)
tn, fp, fn, tp = confusion_matrix(true_class,pred_class).ravel()
print("{}\t{}\t{}\t{}".format(tp,fp,tn,fn))

#END OF PROGRAM
