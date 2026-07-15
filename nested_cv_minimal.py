"""
A simple, leakage-free way to test the DCIS vs Normal classifier.
 
WHAT PROBLEM THIS SOLVES (read this first)
--------------------------------------------
Your original code picked the "best" 1000 CpGs using ALL 55 samples, THEN
split into train/test. That's a problem because picking those CpGs already
looked at the labels of samples that were later used to test the model --
like a student picking exam questions after seeing the answer key.
 
THE FIX
-------
Put the CpG-picking step and the model INSIDE the same "Pipeline" object.
Then, when we cross-validate that Pipeline, the CpG-picking step only ever
gets to see the training data for whichever fold is currently running. It
never sees the test fold's labels. That's the whole fix.
 
HOW TO RUN THIS
----------------
1. Make sure this file and your CSV are in the same folder.
2. Change CSV_PATH below if your file has a different name.
3. Run:  python nested_cv_minimal.py
"""
 
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score
 
# ---------------------------------------------------------------
# STEP 0: Settings you might want to change
# ---------------------------------------------------------------
CSV_PATH = "beta_variance_prefiltered.csv"   # your input file
LABEL_COL = "label"                          # 1 = DCIS, 0 = Normal
NUM_CPGS_TO_KEEP = 500                       # how many CpGs the model gets to use
NUM_FOLDS = 5                                # how many train/test splits to average over
 
# ---------------------------------------------------------------
# STEP 1: Load the data
# ---------------------------------------------------------------
df = pd.read_csv(CSV_PATH, index_col=0)
y = df[LABEL_COL].astype(int).values   # the answers (1 = DCIS, 0 = Normal)
X = df.drop(columns=[LABEL_COL])       # everything else (the CpG beta values)
print(f"Loaded {X.shape[0]} samples and {X.shape[1]} CpGs")
 
# ---------------------------------------------------------------
# STEP 2: Build the pipeline
#   This bundles two steps into one object:
#     (a) SelectKBest -- picks the CpGs that best separate the two groups
#     (b) RandomForestClassifier -- the actual model
#   Because they're bundled, whenever we fit this pipeline on some data,
#   step (a) can ONLY see that same data -- never anything held out.
# ---------------------------------------------------------------
pipeline = Pipeline([
    ("pick_best_cpgs", SelectKBest(score_func=f_classif, k=NUM_CPGS_TO_KEEP)),
    ("model", RandomForestClassifier(
        n_estimators=300,
        class_weight="balanced",   # accounts for 40 DCIS vs 15 Normal
        random_state=42,
    )),
])
 
# ---------------------------------------------------------------
# STEP 3: Test it fairly using 5-fold cross-validation
#   We split the 55 samples into 5 groups ("folds"). One fold at a time
#   becomes the test set, the other 4 become training data. We repeat
#   this 5 times so every sample gets tested exactly once.
# ---------------------------------------------------------------
folds = StratifiedKFold(n_splits=NUM_FOLDS, shuffle=True, random_state=42)
 
fold_scores = []
fold_number = 1
 
for train_rows, test_rows in folds.split(X, y):
    X_train, X_test = X.iloc[train_rows], X.iloc[test_rows]
    y_train, y_test = y[train_rows], y[test_rows]
 
    # Fit the WHOLE pipeline (CpG selection + model) on training data only
    pipeline.fit(X_train, y_train)
 
    # Get predicted probabilities for the held-out test fold
    predicted_probabilities = pipeline.predict_proba(X_test)[:, 1]
 
    # Score this fold
    fold_auc = roc_auc_score(y_test, predicted_probabilities)
    fold_scores.append(fold_auc)
    print(f"Fold {fold_number}: AUC = {fold_auc:.3f}")
    fold_number += 1
 
# ---------------------------------------------------------------
# STEP 4: Report the honest, leakage-free result
# ---------------------------------------------------------------
average_auc = sum(fold_scores) / len(fold_scores)
print(f"\nAverage AUC across all 5 folds: {average_auc:.3f}")
print("This number is safe to report -- CpG selection never saw test data.")
 
# ---------------------------------------------------------------
# STEP 5 (optional): Which CpGs mattered most?
#   We fit the pipeline one final time on ALL the data, just to look at
#   which CpGs it picked and how important each one was. This is only
#   for understanding the model -- NOT a performance score.
# ---------------------------------------------------------------
pipeline.fit(X, y)
chosen_cpgs = X.columns[pipeline.named_steps["pick_best_cpgs"].get_support()]
importance_scores = pipeline.named_steps["model"].feature_importances_
 
top_cpgs = pd.Series(importance_scores, index=chosen_cpgs).sort_values(ascending=False)
print("\nTop 10 most important CpGs:")
print(top_cpgs.head(10))
 
# Save the full list of chosen CpGs + their importance scores to a CSV file
top_cpgs.to_csv("cpg_importance.csv", header=["importance"])
print("\nSaved full CpG importance list to cpg_importance.csv")
