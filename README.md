# 🧬 DNA Methylation-Based ML Classifier
### Distinguishing DCIS from Adjacent-Normal Breast Tissue

`Python` `R` `scikit-learn` `Dataset: GSE66313` `AUC: 0.94`

---

## 📌 Table of Contents
- [Why I Built This](#why-i-built-this)
- [Dataset](#the-dataset)
- [Pipeline Overview](#pipeline-overview)
- [Methods](#methods)
- [Results](#results)
- [Key Findings](#key-findings)
- [A Note on a Bug I Found and Fixed](#a-note-on-a-bug-i-found-and-fixed)
- [Limitations](#limitations)
- [How to Reproduce](#reproducing-this-analysis)
- [Where I'd Take This Next](#where-id-take-this-next)

---

## Why I Built This

This project is a direct extension of my DNA methylation differential analysis work — [DNA-Methylation-Differential-Analysis](#). That project identifies significantly differentially methylated CpG sites between DCIS and Adjacent-Normal breast tissue using R/Bioconductor. This repo takes it one step further: using those CpG sites as features to train a machine learning classifier in Python.

The question I wanted to answer was: can the methylation differences identified statistically in R actually be used to build a model that classifies tissue type? Once the baseline model worked, I pushed it further — including catching and fixing a data leakage bug that was quietly inflating my results (more on that below).

## The Dataset

| Property | Detail |
|---|---|
| Source | GSE66313 — NCBI GEO |
| Platform | Illumina HumanMethylation450 BeadChip (~485,000 CpG sites) |
| Total Samples | 55 breast tissue samples |
| DCIS | 40 samples (label = 1) |
| Adjacent-Normal | 15 samples (label = 0) |

GSE66313 is small enough that every modelling decision is visible and interpretable, but real enough that the results are meaningful — though its small size (55 samples) is exactly what made the leakage bug below so easy to miss and so important to fix.

## Pipeline Overview

This project runs in two stages — a statistical/preprocessing layer in R, then a machine learning layer in Python.

```
Raw IDAT files (Illumina 450K)
          ↓
  ┌─────────────────────────────┐
  │  R — minfi + limma           │
  │  • Background correction     │
  │  • Normalisation              │
  │  • Quality filtering          │
  │  • Differential methylation   │
  │    analysis (for biology,     │
  │    not for feature selection) │
  │  • Export variance-ranked     │
  │    CpGs (label-blind)         │
  └─────────────────────────────┘
          ↓  beta_variance_prefiltered.csv
  ┌──────────────────────────────┐
  │  Python — scikit-learn       │
  │  • CpG selection happens      │
  │    INSIDE cross-validation    │
  │  • Random Forest Classifier   │
  │  • 5-fold cross-validation    │
  │  • Feature importance         │
  └──────────────────────────────┘
```

The file `beta_variance_prefiltered.csv` is the bridge between R and Python — a 55 × 20,001 matrix of beta values for the top 20,000 most variable CpG sites, plus a label column. Either layer of this pipeline can be run independently.

**Why 20,000 variance-ranked CpGs instead of 1,000 differentially-methylated ones?** This is the fix for the leakage bug — explained fully in [its own section below](#a-note-on-a-bug-i-found-and-fixed). Short version: variance doesn't use the DCIS/Normal label at all, so it can't leak it. The label-aware narrowing down to the CpGs that actually matter now happens safely inside Python instead.

## Methods

### Step 1 — Preprocessing and Differential Methylation Analysis (R)

Before touching any machine learning, I ran the standard methylation array QC and normalisation pipeline:

- **minfi** — background correction, quantile normalisation, detection p-value filtering, sex chromosome removal, cross-reactive and polymorphic probe removal (485,512 probes → 351,166 after all filters)
- **limma** — linear model fitted at each CpG site, testing for differential methylation between DCIS and Normal
- **Benjamini-Hochberg FDR correction** for multiple testing across all 351,166 probes

This differential analysis (DMPs, DMRs, GO enrichment) is still a fully valid and useful part of the project — it answers a *biological* question ("which CpGs and pathways differ between DCIS and Normal, using all the data I have?"). It's a completely different question from the *machine learning* one ("can I build a model that generalises to a new, unseen sample?"), and mixing the two up is exactly where the leakage bug crept in.

### Step 2 — Building an Honest Feature Matrix for ML (R)

Instead of picking the top 1000 CpGs by their limma p-value (which uses every sample's label), I rank all 351,166 probes by **variance** across the 55 samples and keep the top 20,000. Variance never looks at which samples are DCIS vs Normal, so this step can't leak label information into the machine learning stage — it's just "which probes vary enough across people to be worth considering at all."

### Step 3 — Machine Learning Classifier (Python)

**The core idea:** CpG selection and the classifier are bundled into one scikit-learn `Pipeline`, so every time the pipeline is fit on a fold of data, the CpG-selection step can only ever see that fold's own training labels — never the labels of samples being tested.

```python
Pipeline([
    ("pick_best_cpgs", SelectKBest(score_func=f_classif, k=500)),
    ("model", RandomForestClassifier(
        n_estimators=300,
        class_weight="balanced",   # accounts for 40 DCIS vs 15 Normal
        random_state=42,
    )),
])
```

**Why Random Forest?**
- Handles high-dimensional data (thousands of candidate CpGs, 55 samples) well
- Produces interpretable feature importance scores
- Natively supports class weighting — critical here, since without `class_weight='balanced'` the model could hit 72.7% accuracy by always predicting DCIS without learning anything about Normal tissue

**Why 5-fold cross-validation, and why it's the whole evaluation this time (no separate train/test split):** with only 55 samples, a single train/test split means the test set is tiny (11 samples — one wrong prediction moves accuracy by 9%). 5-fold CV tests every sample exactly once across five independent rounds, which is both a more honest and a more data-efficient way to evaluate a small dataset. It also means there's no separate "Step 2 model" to accidentally leak into — the whole evaluation lives inside one clean loop.

## Results

### Performance Summary

| Metric | Value |
|---|---|
| Fold 1 AUC | 0.917 |
| Fold 2 AUC | 0.958 |
| Fold 3 AUC | 1.000 |
| Fold 4 AUC | 0.875 |
| Fold 5 AUC | 0.958 |
| **5-Fold CV AUC (mean ± std)** | **0.942 ± 0.042** |

This is the leakage-free number, and the one I'd report going forward. It's a meaningful drop from the 0.983 I originally reported — see [the note below](#a-note-on-a-bug-i-found-and-fixed) for why, and why I think this lower number is actually the more trustworthy result.

Fold 4's dip to 0.875 (against a 1.000 on Fold 3) is a reminder of how much a single sample can move the score with only ~11 samples per test fold. That spread is exactly why averaging across 5 folds, rather than trusting any one split, matters here.

### Feature Importance — Top 10 CpG Sites

| CpG | Importance |
|---|---|
| cg06898823 | 0.0234 |
| cg24280925 | 0.0215 |
| cg01298514 | 0.0199 |
| cg24833737 | 0.0192 |
| cg02006107 | 0.0189 |
| cg26591066 | 0.0187 |
| cg01369207 | 0.0160 |
| cg17928286 | 0.0131 |
| cg02710296 | 0.0129 |
| cg20995977 | 0.0125 |

Encouragingly, several of these sites — `cg06898823`, `cg02006107`, `cg01369207`, `cg17928286` — also came out on top in a separately-tuned run of the same pipeline. `cg17928286` also happens to be the single most significant CpG in the original limma differential methylation analysis. Independent runs and methods landing on overlapping CpGs is a much stronger signal than any one of them alone.

## Key Findings

- **CV AUC 0.942 ± 0.042** — strong, honestly-measured separation of DCIS from Adjacent-Normal tissue
- **No single train/test split is trustworthy at this sample size** — Fold 4 alone would have suggested a much weaker model (AUC 0.875) than Fold 3 (AUC 1.000); only the average across folds is meaningful
- **A handful of CpGs consistently matter across runs** — `cg06898823`, `cg02006107`, `cg01369207`, and `cg17928286` show up near the top of Random Forest importance both in this simplified script and in a separately-tuned run, and `cg17928286` also tops the original limma differential test
- **500 CpGs (out of 20,000 candidates) was enough** — the model doesn't need a huge feature set to perform well, which is encouraging for eventually building a smaller, targeted clinical panel

## A Note on a Bug I Found and Fixed

My original version of this project picked the top 1000 CpGs using a limma test computed on **all 55 samples**, and only split into train/test *after* that. That's data leakage: the feature-selection step already "saw" the labels of samples that later ended up being used to test the model — like a student picking their own exam questions after seeing the answer key.

This inflated my original reported CV AUC to 0.983. It wasn't a wrong number exactly — the model really was finding strong signal — but it wasn't an honest estimate of how well the model would do on a genuinely new sample, because the CpG list itself had already been chosen with knowledge of who was DCIS and who was Normal.

**The fix:** move CpG selection *inside* the cross-validation loop, using an sklearn `Pipeline`, so it only ever sees training-fold labels. I also switched the upstream R export from a label-aware DMP ranking to a label-blind variance ranking, so no step anywhere in the pipeline touches the label before cross-validation does.

The corrected number is **0.942 ± 0.042** — lower than 0.983, and that drop is the size of the leakage I was carrying before. I'm including this section because I think how a bug gets found and fixed is as useful to show as the result itself, especially for anyone learning to build these pipelines themselves.

## Limitations

- **Small N (55 samples):** Cross-validation helps make the most of a small dataset, but this needs validation in a larger, independent cohort before any clinical interpretation.
- **Bulk tissue methylation:** Beta values average across all cell types — epithelial, stromal, immune. Some signal could reflect cell composition differences rather than true tumour methylation. Cell-type deconvolution (e.g., EpiDISH) was not applied.
- **Feature importance method:** This version uses Random Forest's built-in importance (MDI), which is simple to compute but known to be biased for continuous features like methylation beta values. A more rigorous importance ranking would need a different method entirely.
- **No hyperparameter tuning in the simplified script:** `k=500` CpGs and `n_estimators=300` are fixed values rather than grid-searched. A separate, tuned version was tested and gave a nearly identical result (0.950 ± 0.049 AUC), suggesting the model isn't especially sensitive to these choices at this sample size — but it wasn't formally re-verified here.

## Reproducing This Analysis

**R packages:**
```r
BiocManager::install(c("minfi", "limma", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                        "IlluminaHumanMethylation450kmanifest", "missMethyl", "DMRcate",
                        "maxprobes"))
install.packages(c("RColorBrewer", "stringr", "GEOquery", "ggplot2", "EnhancedVolcano", "knitr"))
```

**Python packages:**
```
pip install pandas scikit-learn
```

**Run:**
```r
# Step 1 — R: run dna_methylation_analysis.R through the DMP/DMR/GO analysis,
# then run the final export block to generate the ML feature matrix:
#   results/beta_variance_prefiltered.csv
```
```
# Step 2 — Python: cross-validate the classifier
python nested_cv_minimal.py
```
The Python script expects `beta_variance_prefiltered.csv` in the same folder (edit `CSV_PATH` at the top of the script if yours is named differently). It prints each fold's AUC as it runs, then the final averaged score, then saves the full CpG importance list to `cpg_importance.csv`.

## Project Structure

```
├── dna_methylation_analysis.R      # R pipeline: preprocessing + DMP/DMR/GO + ML export
├── nested_cv_minimal.py            # Python: leakage-free 5-fold CV classifier
├── beta_variance_prefiltered.csv   # Feature matrix — R to Python bridge (leak-free)
├── cpg_importance.csv              # Output: all selected CpGs ranked by importance
├── DMPs_DCIS_vs_Normal.csv         # Output: full differential methylation results
├── DMRs_DCIS_vs_Normal.csv         # Output: differentially methylated regions
├── GO_enrichment_DCIS_vs_Normal.csv# Output: gene ontology enrichment
└── README.md
```

## Where I'd Take This Next

The immediate next step is validation on an independent DCIS cohort — the model performs well on GSE66313, but that needs to hold up on data it has never seen. After that, I'd want to functionally annotate the top CpG sites — map them to gene bodies and promoters, find out which genes sit near them, and check whether those genes have known roles in breast cancer biology.

The longer-term direction that genuinely interests me is testing whether these same CpG sites carry signal in cell-free DNA from patient plasma. That would mean a completely different experimental approach — collecting plasma from DCIS patients, extracting cfDNA, and running bisulfite sequencing targeting these specific loci rather than using an array. Early cancer detection from a blood sample is where this kind of methylation research ultimately points.

---
*Dataset: GSE66313 | Platform: Illumina HumanMethylation450 | Tools: R (minfi, limma) + Python (scikit-learn) | Model: Random Forest*
