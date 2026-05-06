**DNA Methylation-Based Machine Learning Classifier**
Distinguishing DCIS from Adjacent-Normal Breast Tissue | GSE66313 | Illumina 450K

**Why I Built This**
This project is a direct extension of my DNA methylation differential analysis work — DNA-Methylation-Differential-Analysis. That project identifies significantly differentially methylated CpG sites between DCIS and Adjacent-Normal breast tissue using R/Bioconductor. This repo takes it one step further: using those top CpG sites as features to train a machine learning classifier in Python.
The question I wanted to answer was: can the methylation differences identified statistically in R actually be used to build a model that classifies tissue type? And once the baseline model worked, I wanted to push it — run hyperparameter tuning, compare models, and understand which specific CpG sites were driving the classification.

# The Dataset
Source: GSE66313 — publicly available on NCBI GEO
Platform: Illumina HumanMethylation450 BeadChip (~485,000 CpG sites)
Samples: 55 breast tissue samples

# 40 DCIS (Ductal Carcinoma In Situ) — label = 1
15 Adjacent-Normal — label = 0

The samples are paired — each patient contributed both tumour and adjacent-normal tissue. That design matters: it means the methylation differences I found are more likely to reflect true cancer biology rather than normal variation between people.
GSE66313 is small enough that every modelling decision is visible and interpretable, but real enough that the results are actually meaningful. It's not a toy dataset.

# Pipeline Overview
This project runs in two stages — a statistical layer in R to find signal, then a machine learning layer in Python to build a classifier on that signal. That two-step structure reflects how this kind of analysis actually works in practice.
Raw IDAT files (Illumina 450K)
        ↓
  R — minfi + limma
  • Background correction (ssNoob)
  • Normalisation (quantile)
  • Quality filtering
  • Differential methylation analysis
  • Export top 1000 DMPs as beta matrix
        ↓
  Python — scikit-learn
  • Stratified train/test split
  • Random Forest Classifier
  • 5-fold stratified cross-validation
  • Hyperparameter tuning (GridSearchCV)
  • Feature importance analysis
  • Visualisation
The file connecting R and Python is beta_top1000_DMPs.csv — a 55 × 1001 matrix of beta values for the top 1000 differentially methylated CpG sites, plus a label column. Either layer of this pipeline can be run independently.

# Methods
# Step 1 — Differential Methylation Analysis (R)
Before touching any machine learning, I first needed to figure out which CpG sites actually differ between DCIS and Normal tissue. Running all 485,000 CpGs through a classifier on 55 samples would be a recipe for overfitting — you'd have nearly 9,000 features per sample, which is a statistical disaster.
Instead, I used minfi for preprocessing (background correction, normalisation, quality filtering) and limma for differential methylation testing — fitting a linear model at each CpG site and adjusting for multiple testing using Benjamini-Hochberg FDR correction. The top 1000 sites ranked by adjusted p-value became the feature set for the classifier. These aren't random features — they're the sites with the strongest statistical evidence of genuine methylation differences between the two tissue types.
# Step 2 — Machine Learning Classifier (Python)
Feature matrix: 55 samples × 1000 CpG beta values. Each beta value is between 0 and 1, representing the proportion of methylated alleles at that site across the sample.
Train/test split: 80/20 stratified split — 44 samples for training, 11 for testing. Stratification makes sure both splits maintain the 40:15 DCIS:Normal ratio, so the model doesn't end up training on only one class by bad luck.
Model:
pythonRandomForestClassifier(
    n_estimators=500,        # 500 trees voting by majority
    max_features='sqrt',     # each tree sees ~32 CpGs — keeps trees diverse
    class_weight='balanced', # corrects for 40 DCIS vs 15 Normal imbalance
    random_state=42,
    n_jobs=-1
)
I chose Random Forest because it handles high-dimensional data well, gives feature importance scores, and lets you correct for class imbalance directly. With 40 DCIS and only 15 Normal samples, that last point matters — without class_weight='balanced', the model could hit 72.7% accuracy just by always predicting DCIS, without learning anything useful about Normal tissue.
Why cross-validation over a single test split? The test set is only 11 samples. One wrong prediction equals a 9% accuracy change — that's not a stable number to report. Five-fold stratified cross-validation tests every sample exactly once across five independent rounds, giving a much more honest picture of how the model actually performs.
Hyperparameter tuning: After getting the baseline results, I ran GridSearchCV over 18 parameter combinations (each evaluated with 5-fold CV) to check whether the baseline could be improved.

# Results
MetricValueTest Accuracy81.8% (9/11 correct)Test AUC-ROC0.9175-Fold CV AUC0.983 ± 0.0205-Fold CV Accuracy94.5% ± 4.5%
The CV AUC of 0.983 with a standard deviation of ±0.020 tells me the model is consistently separating DCIS from Normal across all five folds — it's not a one-fold fluke. The test accuracy of 81.8% looks lower, but that's a small-N artefact. With 11 test samples, the cross-validated 94.5% is the number I trust.
Baseline vs GridSearch — Did Tuning Help?
When I got CV AUC 0.983 with the baseline model, my first instinct was to try to push it further. GridSearchCV tested all 18 combinations of these parameters, each with 5-fold CV — 90 models in total:
pythonparams = {
    'n_estimators':      [100, 300, 500],
    'max_depth':         [None, 5, 10],
    'min_samples_split': [2, 5]
}
Baseline RFTuned RF (GridSearch)n_estimators500300max_depthNoneNonemin_samples_split22Test Accuracy0.8180.818Test AUC0.9170.917CV AUC (mean ± std)0.983 ± 0.0200.983 ± 0.020CV Accuracy (mean ± std)0.945 ± 0.0450.945 ± 0.045
The tuned model matched the baseline exactly. Honestly, I was hoping for an improvement — but what I got instead was something more useful: confirmation that the baseline wasn't a lucky guess. The settings were already right for this dataset. The one interesting detail is that 300 trees performed identically to 500, meaning the model stabilises early. You don't need a heavier ensemble to get the same result.
The fact that GridSearch couldn't beat the baseline also means the 0.983 AUC isn't sensitive to parameter choices — which reduces the worry that it's an inflated number.
# Confusion Matrix
Show Image
9 out of 11 correct. The model missed 1 Normal (called it DCIS) and 1 DCIS (called it Normal). The missed DCIS is the more clinically concerning error — a false negative means a missed cancer. But with only 3 Normal samples in the test set, it's hard to draw strong conclusions from single errors at this scale.
# ROC Curve
Show Image
The curve sits close to the top-left corner, which is where you want it. The stepped shape is just an artefact of having few test samples — not a sign of instability. The 5-fold CV AUC of 0.983 is the more reliable number.
Feature Importance — Top 20 CpG Sites
Show Image
The top 5 CpGs carry noticeably more weight than the remaining 15. That concentration matters — it suggests the signal separating DCIS from Normal isn't spread evenly across all 1000 sites. A small targeted panel of CpGs might be enough to build a viable classifier, which is relevant if you're thinking about translating this into a clinical assay.
Methylation Distribution — Top 4 CpGs
Show Image

cg17302155 — higher in DCIS (~0.25) than Normal (~0.10). Hypermethylated in cancer, consistent with silencing of a tumour suppressor gene.
cg03550233 — lower in DCIS (~0.18) than Normal (~0.26). Hypomethylated, which can activate oncogenes or destabilise repeat elements.
cg14717069 — lower in DCIS (~0.26) than Normal (~0.40). Another hypomethylation pattern in cancer tissue.
cg13643356 — higher in DCIS (~0.57) than Normal (~0.40), but with a wide spread in DCIS — suggesting this site varies quite a bit between individual tumours.

The fact that some CpGs go up and some go down in DCIS is biologically expected. Cancer methylation reprogramming is site-specific, not a global shift in one direction. Seeing both patterns in the top features adds confidence that the classifier is picking up real signal.

# Key Findings
The model achieves CV AUC 0.983 — near-perfect separation of DCIS from Adjacent-Normal tissue using 1000 methylation features. The signal isn't spread evenly across those 1000 sites; it's concentrated in a handful of CpGs, with the top 5 carrying disproportionate weight. GridSearchCV confirmed the baseline was already optimal — 90 models, no improvement. And the direction of methylation change at the top CpGs (some up, some down in DCIS) is consistent with what we know about cancer epigenetics, which gives the results biological credibility beyond just the numbers.

# Limitations
No analysis is perfect, and I'd rather flag these upfront than have someone else point them out:

Small N (55 samples): The test set is only 11 samples. Cross-validation helps, but this needs validation in a larger independent cohort before any clinical interpretation.
Bulk tissue methylation: Beta values average across all cell types in the sample — epithelial, stromal, immune. Some of the signal could reflect differences in cell composition between DCIS and Normal rather than true tumour methylation changes. Cell-type deconvolution (e.g., EpiDISH) was not applied.
Feature selection on the full dataset: I selected the top 1000 DMPs using all 55 samples before splitting. Strictly, this should happen inside the training fold to avoid any information leakage. It's a known limitation I'd fix in a follow-up.
Feature importance method: MDI from Random Forest can overestimate importance for continuous features. SHAP values would give a more reliable ranking.


# Reproducing This Analysis
R packages:
rBiocManager::install(c("minfi", "limma", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
Python packages:
bashpip install pandas scikit-learn matplotlib seaborn
Run:
bash# Step 1 — R: generates beta_top1000_DMPs.csv
Rscript dna_methylation_analysis.R

# Step 2 — Python: generates ml_results/ with all plots and metrics
python DNA_methylation_ML_Classifier.py

Project Structure
├── dna_methylation_analysis.R           # R pipeline: preprocessing + differential methylation
├── DNA_methylation_ML_Classifier.py     # Python ML classifier + GridSearchCV + plots
├── beta_top1000_DMPs.csv                # Feature matrix exported from R — bridge to Python
├── ml_results/
│   ├── roc_curve.png
│   ├── confusion_matrix.png
│   ├── feature_importance.png
│   ├── top4_cpg_boxplots.png
│   ├── top20_features.csv
│   └── model_metrics.csv
└── README.md

Note: The R script (dna_methylation_analysis.R) is the same pipeline from the parent project DNA-Methylation-Differential-Analysis. Running it on GSE66313 IDAT files produces beta_top1000_DMPs.csv, which is the input to the Python classifier.


**Where I'd Take This Next**
The immediate next step is validation on an independent DCIS cohort — the model performs well on GSE66313, but that needs to hold up on data it has never seen. After that, I'd want to functionally annotate the top CpG sites — map them to gene bodies and promoters, find out which genes sit near them, and check whether those genes have known roles in breast cancer biology.
The longer-term direction that genuinely interests me is testing whether these same CpG sites carry signal in cell-free DNA from patient plasma. That would mean a completely different experimental approach — collecting plasma from DCIS patients, extracting cfDNA, and running bisulfite sequencing targeting these specific loci rather than using an array. The ML classification step would follow from that. It's a full new study rather than a direct extension, but it's the direction that makes this work clinically relevant. Early cancer detection from a blood sample is where this kind of methylation research ultimately points.

Dataset: GSE66313 | Platform: Illumina HumanMethylation450 | Tools: R (minfi, limma) + Python (scikit-learn) | Model: Random Forest
