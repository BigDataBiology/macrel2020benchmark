import pandas as pd
from sklearn import metrics, ensemble
from os import makedirs
import pickle
import gzip

makedirs('models/', exist_ok=True)
def report(name, y_true, y_pred):
    print(f"# {name} CLASSIFIER ")
    print("Accuracy: {:.3}".format(metrics.accuracy_score(y_true, y_pred)))
    print("MCC: {:.3}".format(metrics.matthews_corrcoef(y_true, y_pred)))

    print()
    print("Confusion matrix:")
    print(metrics.confusion_matrix(y_true, y_pred))
    print()
    print(metrics.classification_report(y_true, y_pred, digits=3))
    print("\n")

train11 = pd.read_table('1_1_trainset.tsv', index_col=0)
train15 = pd.read_table('1_5_trainset.tsv', index_col=0)
train110 = pd.read_table('1_10_trainset.tsv', index_col=0)
train120 = pd.read_table('1_20_trainset.tsv', index_col=0)
train130 = pd.read_table('1_30_trainset.tsv', index_col=0)
train140 = pd.read_table('1_40_trainset.tsv', index_col=0)
train150 = pd.read_table('1_50_trainset.tsv', index_col=0)
test11 = pd.read_table('test_set_1_1.tsv', index_col=0)


rf_w_oob = ensemble.RandomForestClassifier(
        n_estimators=101,
        oob_score=True,
        random_state=12345,
        n_jobs=8)

# We train a model without the OOB predictions to save to disk
rf_wout_oob = ensemble.RandomForestClassifier(
        n_estimators=101,
        random_state=12345,
        n_jobs=8)
rf_w_oob.fit(train150.iloc[:, 3:], train150['group'])
rf_wout_oob.fit(train150.iloc[:, 3:], train150['group'])

# In testing, using joblib.dump was actually slightly worse than standard pickle in
# terms of diskspace
pickle.dump(rf_wout_oob, gzip.open('models/AMP.pkl.gz', 'wb'))


oob_pred = rf_w_oob.classes_[(rf_w_oob.oob_decision_function_.T[1] > .5).astype(int)]
test_pred = rf_w_oob.predict(test11.iloc[:, 3:])

# As xiao_test & ampep_test overlap in sequences, we after clustering sequences
# do not need to do out-of-bag predictions:
rf_wout_oob.fit(train11.iloc[:,3:], train11['group'])
test_pred = rf_wout_oob.predict(test11.iloc[:, 3:])
report('(1:1 AMP:NAMP)', test11['group'], test_pred)

rf_wout_oob.fit(train15.iloc[:,3:], train15['group'])
test_pred = rf_wout_oob.predict(test11.iloc[:, 3:])
report('(1:5 AMP:NAMP)', test11['group'], test_pred)

rf_wout_oob.fit(train110.iloc[:,3:], train110['group'])
test_pred = rf_wout_oob.predict(test11.iloc[:, 3:])
report('(1:10 AMP:NAMP)', test11['group'], test_pred)

rf_wout_oob.fit(train120.iloc[:,3:], train120['group'])
test_pred = rf_wout_oob.predict(test11.iloc[:, 3:])
report('(1:20 AMP:NAMP)', test11['group'], test_pred)

rf_wout_oob.fit(train130.iloc[:,3:], train130['group'])
test_pred = rf_wout_oob.predict(test11.iloc[:, 3:])
report('(1:30 AMP:NAMP)', test11['group'], test_pred)

rf_wout_oob.fit(train140.iloc[:,3:], train140['group'])
test_pred = rf_wout_oob.predict(test11.iloc[:, 3:])
report('(1:40 AMP:NAMP)', test11['group'], test_pred)

rf_wout_oob.fit(train150.iloc[:,3:], train150['group'])
test_pred = rf_wout_oob.predict(test11.iloc[:, 3:])
report('(1:50 AMP:NAMP)', test11['group'], test_pred)
