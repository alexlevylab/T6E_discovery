import glob 

import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, confusion_matrix
from sklearn.model_selection import cross_val_score
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score

np.random.seed(613)

#######

dfs = []

for file in glob.glob("../train_data/scores_maxima.csv"):
	dfs.append(pd.read_csv(file))

df = pd.concat(dfs)

#######
scaler = StandardScaler()
X = scaler.fit_transform(df['mean_intermolecular_PAE'].to_numpy().reshape(-1, 1))
y = df['expected'].to_numpy().reshape(-1, 1).ravel()

X_train, X_test, y_train, y_test = train_test_split(X,y, test_size = 0.3)


#####
print("\nlogistic regression")

model = LogisticRegression(C=0.25)

fivefold_cv = cross_val_score(model, X_train, y_train, cv = 8, scoring = 'accuracy')
print(fivefold_cv)

model.fit(X_train, y_train)
pred = model.predict(X_test)

cm = confusion_matrix(y_test, pred)
acc = accuracy_score(y_test, pred)
f1 = f1_score(y_test, pred)
pre = precision_score(y_test, pred)
rec = recall_score(y_test, pred)

print(cm)
print(acc)
print(f1)
print(pre)
print(rec)

fpr, tpr, thresholds = roc_curve(y_test, pred)
AUC = roc_auc_score(y_test, pred)

##
prob_true = model.predict_proba(X)[:,1]
X_unscaled = scaler.inverse_transform(X)

#####

print("\nreal data from paper")

test = pd.read_csv("../test_data/scores_maxima.csv")

scaler = StandardScaler()
X_validation = scaler.fit_transform(test['mean_intermolecular_PAE'].to_numpy().reshape(-1, 1))
y_test = test['expected'].to_numpy().reshape(-1, 1)

pred = model.predict(X_validation)

cm = confusion_matrix(y_test, pred)
acc = accuracy_score(y_test, pred)
f1 = f1_score(y_test, pred)
pre = precision_score(y_test, pred)
rec = recall_score(y_test, pred)

print(cm)
print(acc)
print(f1)
print(pre)
print(rec)



