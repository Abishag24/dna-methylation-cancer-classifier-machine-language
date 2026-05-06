# 🧬 DNA Methylation-Based ML Classifier
### Distinguishing DCIS from Adjacent-Normal Breast Tissue

![Python](https://img.shields.io/badge/Python-3.8+-blue?logo=python)
![R](https://img.shields.io/badge/R-Bioconductor-276DC3?logo=r)
![scikit-learn](https://img.shields.io/badge/scikit--learn-RandomForest-orange)
![Dataset](https://img.shields.io/badge/Dataset-GSE66313-green)
![AUC](https://img.shields.io/badge/CV%20AUC-0.983-brightgreen)

---

## 📌 Table of Contents
- [Why I Built This](#why-i-built-this)
- [Dataset](#the-dataset)
- [Pipeline Overview](#pipeline-overview)
- [Methods](#methods)
- [Results](#results)
- [Key Findings](#key-findings)
- [Limitations](#limitations)
- [How to Reproduce](#reproducing-this-analysis)
- [Where I'd Take This Next](#where-id-take-this-next)

---

## Why I Built This

This project is a direct extension of my DNA methylation differential analysis work — [DNA-Methylation-Differential-Analysis](https://github.com/Abishag24/DNA-Methylation-Differential-Analysis). That project identifies significantly differentially methylated CpG sites between DCIS and Adjacent-Normal breast tissue using R/Bioconductor. This repo takes it one step further: using those top CpG sites as features to train a machine learning classifier in Python.

The question I wanted to answer was: **can the methylation differences identified statistically in R actually be used to build a model that classifies tissue type?** And once the baseline model worked, I wanted to push it — run hyperparameter tuning, compare models, and understand which specific CpG sites were driving the classification.

---

## The Dataset

| Property | Detail |
|----------|--------|
| **Source** | [GSE66313](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66313) — NCBI GEO |
| **Platform** | Illumina HumanMethylation450 BeadChip (~485,000 CpG sites) |
| **Total Samples** | 55 breast tissue samples |
| **DCIS** | 40 samples (label = 1) |
| **Adjacent-Normal** | 15 samples (label = 0) |
| **Design** | Paired — each patient contributed both tumour and adjacent-normal tissue |

The paired design matters: it means the methylation differences I found are more likely to reflect true cancer biology rather than normal variation between people. GSE66313 is small enough that every modelling decision is visible and interpretable, but real enough that the results are meaningful. It's not a toy dataset.

---

## Pipeline Overview

This project runs in two stages — a statistical layer in R to find signal, then a machine learning layer in Python to build a classifier on that signal.

```
Raw IDAT files (Illumina 450K)
          ↓
  ┌─────────────────────────┐
  │  R — minfi + limma      │
  │  • Background correction│
  │  • Normalisation        │
  │  • Quality filtering    │
  │  • Differential methyl. │
  │  • Export top 1000 DMPs │
  └─────────────────────────┘
          ↓  beta_top1000_DMPs.csv
  ┌──────────────────────────────┐
  │  Python — scikit-learn       │
  │  • Stratified train/test     │
  │  • Random Forest Classifier  │
  │  • 5-fold cross-validation   │
  │  • GridSearchCV tuning       │
  │  • Feature importance        │
  │  • Visualisation             │
  └──────────────────────────────┘
```

> The file `beta_top1000_DMPs.csv` is the bridge between R and Python — a 55 × 1001 matrix of beta values for the top 1000 differentially methylated CpG sites, plus a label column. Either layer of this pipeline can be run independently.

---

## Methods

### Step 1 — Differential Methylation Analysis (R)

Before touching any machine learning, I needed to identify which CpG sites actually differ between DCIS and Normal tissue. Running all 485,000 CpGs through a classifier on 55 samples would be a recipe for overfitting — nearly 9,000 features per sample.

Instead:
- **`minfi`** — preprocessing: background correction (ssNoob), quantile normalisation, quality filtering
- **`limma`** — linear model fitted at each CpG site, testing for differential methylation
- **Benjamini-Hochberg FDR** correction for multiple testing
- **Top 1000 CpGs** ranked by adjusted p-value selected as features

These aren't random features — they're the sites with the strongest statistical evidence of genuine methylation differences between the two tissue types.

---

### Step 2 — Machine Learning Classifier (Python)

**Data setup:**
- Feature matrix: 55 samples × 1000 CpG beta values (0–1 scale)
- 80/20 stratified train/test split → 44 training, 11 test samples
- Stratification preserves the 40:15 DCIS:Normal ratio in both splits

**Model:**

```python
RandomForestClassifier(
    n_estimators=500,        # 500 trees voting by majority
    max_features='sqrt',     # each tree sees ~32 CpGs — keeps trees diverse
    class_weight='balanced', # corrects for 40 DCIS vs 15 Normal imbalance
    random_state=42,
    n_jobs=-1
)
```

**Why Random Forest?**
- Handles high-dimensional data (1000 features, 55 samples) well
- Produces interpretable feature importance scores
- Natively supports class weighting — critical here, since without `class_weight='balanced'`, the model could hit 72.7% accuracy by always predicting DCIS without learning anything about Normal tissue

**Why cross-validation?**
The test set is only 11 samples — one wrong prediction equals a 9% accuracy change. Five-fold stratified CV tests every sample exactly once across five independent rounds, giving a much more honest performance estimate.

**Hyperparameter tuning:**
GridSearchCV tested 18 parameter combinations, each evaluated with 5-fold CV — 90 models total.

---

## Results

### Performance Summary

| Metric | Value |
|--------|-------|
| Test Accuracy | 81.8% (9/11 correct) |
| Test AUC-ROC | 0.917 |
| **5-Fold CV AUC** | **0.983 ± 0.020** |
| **5-Fold CV Accuracy** | **94.5% ± 4.5%** |

The CV AUC of 0.983 ± 0.020 tells me the model is consistently separating DCIS from Normal across all five folds — it's not a one-fold fluke. The test accuracy of 81.8% looks lower, but that's a small-N artefact. With only 11 test samples, the cross-validated 94.5% is the number I trust.

---

### Baseline vs GridSearch — Did Tuning Help?

When the baseline gave CV AUC 0.983, my first instinct was to push it further. GridSearchCV tested:

```python
params = {
    'n_estimators':      [100, 300, 500],
    'max_depth':         [None, 5, 10],
    'min_samples_split': [2, 5]
}
```

| | Baseline RF | Tuned RF (GridSearch) |
|---|---|---|
| n_estimators | 500 | 300 |
| max_depth | None | None |
| min_samples_split | 2 | 2 |
| Test Accuracy | 0.818 | 0.818 |
| Test AUC | 0.917 | 0.917 |
| CV AUC (mean ± std) | 0.983 ± 0.020 | 0.983 ± 0.020 |
| CV Accuracy (mean ± std) | 0.945 ± 0.045 | 0.945 ± 0.045 |

The tuned model matched the baseline exactly. Honestly, I was hoping for an improvement — but what I got instead was confirmation that the baseline wasn't a lucky guess. The settings were already right for this dataset. The one interesting detail is that 300 trees performed identically to 500 — the model stabilises early and doesn't need a heavier ensemble to get the same result.

---

### Confusion Matrix

![Confusion Matrix](ml_results/confusion_matrix.png)

9 out of 11 correct:
- ✅ 2 Normal correctly identified
- ✅ 7 DCIS correctly identified
- ⚠️ 1 Normal wrongly called DCIS (false positive)
- ❌ 1 DCIS wrongly called Normal (false negative — the more clinically concerning error)

---

### ROC Curve

![ROC Curve](ml_results/roc_curve.png)

The curve sits close to the top-left corner — the ideal geometry for a binary classifier. The stepped shape is an artefact of small test sets, not model instability. The 5-fold CV AUC of 0.983 is the more reliable number.

---

### Feature Importance — Top 20 CpG Sites

![Feature Importance](ml_results/feature_importance.png)

The top 5 CpGs carry noticeably more weight than the remaining 15. That concentration matters — it suggests the signal separating DCIS from Normal isn't spread evenly. A small targeted panel of CpGs might be enough to build a viable classifier, which is relevant for translating this into a clinical assay.

---

### Methylation Distribution — Top 4 CpGs

![Methylation Boxplots](ml_results/top4_cpg_boxplots.png)

| CpG Site | Normal | DCIS | Direction | Interpretation |
|----------|--------|------|-----------|----------------|
| cg17302155 | ~0.10 | ~0.25 | ⬆️ Hypermethylated | Consistent with tumour suppressor silencing |
| cg03550233 | ~0.26 | ~0.18 | ⬇️ Hypomethylated | May activate oncogenes or destabilise repeat elements |
| cg14717069 | ~0.40 | ~0.26 | ⬇️ Hypomethylated | Aberrant demethylation in cancer tissue |
| cg13643356 | ~0.40 | ~0.57 | ⬆️ Hypermethylated | Wide DCIS spread suggests inter-tumour heterogeneity |

The mix of hyper- and hypomethylated sites is biologically expected — cancer methylation reprogramming is site-specific, not a global shift in one direction. Seeing both patterns in the top features adds confidence the classifier is picking up real signal.

---

## Key Findings

- **CV AUC 0.983** — near-perfect separation of DCIS from Adjacent-Normal using 1000 methylation features
- **Signal is concentrated** — top 5 CpGs carry disproportionate discriminative weight out of 1000; a smaller panel may be sufficient for a clinical assay
- **GridSearch confirmed baseline optimality** — 90 models, no improvement; the result isn't sensitive to parameter choices
- **Biologically plausible patterns** — top CpGs show both hyper- and hypomethylation in DCIS, consistent with known cancer epigenetic reprogramming

---

## Limitations

No analysis is perfect — I'd rather flag these upfront:

- **Small N (55 samples):** Test set is only 11 samples. Cross-validation helps, but this needs validation in a larger independent cohort before any clinical interpretation.
- **Bulk tissue methylation:** Beta values average across all cell types — epithelial, stromal, immune. Some signal could reflect cell composition differences rather than true tumour methylation. Cell-type deconvolution (e.g., EpiDISH) was not applied.
- **Feature selection on full dataset:** Top 1000 DMPs were selected using all 55 samples before the train/test split — a mild leakage issue I'd fix in a follow-up by moving DMP selection inside the cross-validation loop.
- **Feature importance method:** MDI from Random Forest can overestimate importance for continuous features. SHAP values would give a more reliable ranking.

---

## Reproducing This Analysis

**R packages:**
```r
BiocManager::install(c("minfi", "limma", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
```

**Python packages:**
```bash
pip install pandas scikit-learn matplotlib seaborn
```

**Run:**
```bash
# Step 1 — R: generates beta_top1000_DMPs.csv
Rscript dna_methylation_analysis.R

# Step 2 — Python: generates ml_results/ with all plots and metrics
python DNA_methylation_ML_Classifier.py
```

---

## Project Structure

```
├── dna_methylation_analysis.R           # R pipeline: preprocessing + DMP analysis
├── DNA_methylation_ML_Classifier.py     # Python ML classifier + GridSearchCV + plots
├── beta_top1000_DMPs.csv                # Feature matrix — R to Python bridge
├── ml_results/
│   ├── roc_curve.png
│   ├── confusion_matrix.png
│   ├── feature_importance.png
│   ├── top4_cpg_boxplots.png
│   ├── top20_features.csv
│   └── model_metrics.csv
└── README.md
```

> **Note:** The R script is the same pipeline from the parent project [DNA-Methylation-Differential-Analysis](https://github.com/Abishag24/DNA-Methylation-Differential-Analysis). Running it on GSE66313 IDAT files produces `beta_top1000_DMPs.csv`, which feeds into the Python classifier.

---

## Where I'd Take This Next

The immediate next step is validation on an independent DCIS cohort — the model performs well on GSE66313, but that needs to hold up on data it has never seen. After that, I'd want to functionally annotate the top CpG sites — map them to gene bodies and promoters, find out which genes sit near them, and check whether those genes have known roles in breast cancer biology.

The longer-term direction that genuinely interests me is testing whether these same CpG sites carry signal in cell-free DNA from patient plasma. That would mean a completely different experimental approach — collecting plasma from DCIS patients, extracting cfDNA, and running bisulfite sequencing targeting these specific loci rather than using an array. The ML classification step would follow from that. It's a full new study rather than a direct extension, but it's the direction that makes this work clinically relevant. Early cancer detection from a blood sample is where this kind of methylation research ultimately points.

---

*Dataset: GSE66313 &nbsp;|&nbsp; Platform: Illumina HumanMethylation450 &nbsp;|&nbsp; Tools: R (minfi, limma) + Python (scikit-learn) &nbsp;|&nbsp; Model: Random Forest*
