import pickle
import pandas as pd
from sklearn.metrics import confusion_matrix, accuracy_score, f1_score, precision_score, recall_score, roc_curve, roc_auc_score
import numpy as np

np.random.seed(613)

# Load the trained model from the file
with open('trained_model_iptm_ptm_gteq0.5.pkl', 'rb') as f:
    model = pickle.load(f)

# Load the scaler from the file
with open('scaler_iptm_ptm_gteq0.5.pkl', 'rb') as f:
    scaler = pickle.load(f)

# 'df' with scores and labels
df = pd.read_csv("dataframe.csv")

df = df.loc[df['ptm'] >=0.5]

X_test = df['iptm']

# Replace 'labels_column' with the name of the column containing the labels (ground truth)
y_test = df['expected']

# Apply the scaler to the data
X_test = scaler.transform(X_test.values.reshape(-1, 1))

# Make predictions using the loaded model
pred = model.predict(X_test)
out = pd.concat([df.reset_index(), pd.Series(pred).reset_index()], axis = 1)
out.to_csv("df_withPredictions.csv", index = False)

# Evaluate the model
cm = confusion_matrix(y_test, pred)
acc = accuracy_score(y_test, pred)
f1 = f1_score(y_test, pred)
pre = precision_score(y_test, pred)
rec = recall_score(y_test, pred)

print("Confusion Matrix:")
print(cm)
print("Accuracy:", acc)
print("F1 Score:", f1)
print("Precision:", pre)
print("Recall:", rec)

# Calculate ROC curve and AUC
fpr, tpr, thresholds = roc_curve(y_test, pred)
AUC = roc_auc_score(y_test, pred)

print("AUC:", AUC)

