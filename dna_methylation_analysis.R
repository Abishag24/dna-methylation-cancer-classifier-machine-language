# ---- SECTION 1: LOAD LIBRARIES -----------------------------------------------

#BiocManager::install("maxprobes")


library(knitr)  
library(limma)  
library(minfi)  
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)   # or the annotation manifest matching your array
library(IlluminaHumanMethylation450kmanifest)             # ensure this matches your array version  
library(RColorBrewer)  
library(missMethyl)  
library(DMRcate)  
library(stringr)
library(GEOquery)
library(ggplot2)
library(EnhancedVolcano)
library(maxprobes)


# ---- SECTION 1A: SET WORKING DIRECTORY ----------------------------------------
# Set to your local path where IDAT files are stored
# Modify this path to match your local IDAT file location
setwd("F:/Downloads/GSE66313_RAW")
dataDirectory <- getwd()
cat("Working directory:", dataDirectory, "\n")

options(ExperimentHub.ask = FALSE)                                    # Use current working directory


# ---- SECTION 2: EXTRACT IDAT FILES (run once only, then comment out) ---------

# library(R.utils)
# gz_files <- list.files("G:/Downloads/GSE66313_RAW/",
#                        pattern = "\\.gz$",
#                        full.names = TRUE)
# sapply(gz_files, gunzip, overwrite = TRUE)

# ---- SECTION 3: READ IDAT FILES ----------------------------------------------

rgset <- read.metharray.exp(base = dataDirectory)

# ---- SECTION 4: LOAD AND ALIGN PHENOTYPE DATA FROM GEO ----------------------

gse <- getGEO("GSE66313", GSEMatrix = TRUE, getGPL = FALSE)
pheno <- pData(gse[[1]])

# *** VERIFY YOUR EXACT FACTOR LEVELS BEFORE PROCEEDING ***
cat("\n--- Tissue types in dataset ---\n")
print(table(pheno$`tissue:ch1`))

# Align pheno rows to match the order of samples in rgset
gsm_ids <- sub("_.*", "", sampleNames(rgset))
pheno <- pheno[match(gsm_ids, pheno$geo_accession), ]

# ---- SECTION 5: CREATE GROUP FACTOR ------------------------------------------
# Use tissue:ch1 as the group variable (all 55 samples)
# IMPORTANT: the strings below ("DCIS", "adjacent-normal") must exactly match

group <- factor(pheno$`tissue:ch1`)
levels(group) <- make.names(levels(group))  
pData(rgset)$Group <- group

cat("\n--- Group factor levels (confirm these are correct) ---\n")
print(levels(group))   

# ---- SECTION 6: LOAD 450K ANNOTATION -----------------------------------------
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(anno450k)

# ---- SECTION 7: QUALITY CONTROL — DETECTION P-VALUES ------------------------
detP <- detectionP(rgset)
meanDetP <- colMeans(detP) 
cat("\n--- Mean detection p-values per sample ---\n")
print(meanDetP)

# Identify failed samples (mean detP >= 0.05)
cat("\n--- Samples passing QC (TRUE = pass) ---\n")
print(meanDetP < 0.05)

##---- SECTION 8: QC BOXPLOT --------------------------------------------------

pal <- brewer.pal(8, "Dark2")
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Create a data frame for plotting
detP_df <- data.frame(
  DetP = as.vector(detP),
  Sample = rep(colnames(detP), each = nrow(detP)),
  Group = rep(pData(rgset)$Group, each = nrow(detP))
)

# LEFT PANEL: Boxplot by sample
boxplot(colMeans(detP) ~ pData(rgset)$Group,
        col = pal[1:2],
        ylab = "Mean detection p-value",
        xlab = "Sample Group",
        main = "QC: Mean Detection P-values",
        outline = TRUE,
        las = 1)
abline(h = 0.05, col = "red", lwd = 2, lty = 2)
abline(h = 0.01, col = "blue", lwd = 1, lty = 3)
legend("topright",
       legend = c("p = 0.05", "p = 0.01"),
       lty = c(2, 3),
       col = c("red", "blue"),
       lwd = c(2, 1),
       bg = "white")

# RIGHT PANEL: Individual sample points
plot(colMeans(detP),
     col = pal[factor(pData(rgset)$Group)],
     pch = 19,
     cex = 1.2,
     ylab = "Mean detection p-value",
     xlab = "Sample index",
     main = "Individual Samples",
     ylim = c(0, max(colMeans(detP)) * 1.2))
abline(h = 0.05, col = "red", lwd = 2, lty = 2)
abline(h = 0.01, col = "blue", lwd = 1, lty = 3)
legend("topright",
       legend = levels(factor(pData(rgset)$Group)),
       col = pal[1:2],
       pch = 19,
       bg = "white")

dev.print(png, "QC_DetectionP_Boxplots.png", width = 1200, height = 600, res = 150)
cat("Saved: QC_DetectionP_Boxplots.png\n")

# ---- SECTION 9: REMOVE FAILED SAMPLES ----------------------------------------
keep  <- (meanDetP < 0.05)
rgset <- rgset[, keep]
pheno <- pheno[keep, ]
detP <- detP[, keep]

cat("\nSamples remaining after QC:", ncol(rgset), "\n")

# ---- SECTION 10: NORMALIZATION & QUALITY FILTERING ---------------------------

# 10A: Normalize
mSetSq <- preprocessQuantile(rgset)
mSetRaw <- preprocessRaw(rgset)
cat("\n--- After normalization ---\n")
cat("Probes:", nrow(mSetSq), "\n")

# 10B: Remove failed probes
keep <- rowSums(detP < 0.01) == ncol(detP)
mSetSq <- mSetSq[keep, ]
cat("\n--- After detection p-value filtering ---\n")
cat("Probes remaining:", nrow(mSetSq), "\n")
table(keep)

# 10C: Remove sex chromosomes
sex_cpgs <- anno450k$Name[anno450k$chr %in% c("chrX", "chrY")]
keep <- !(featureNames(mSetSq) %in% sex_cpgs)
mSetSq <- mSetSq[keep, ]
cat("\n--- After sex chromosome removal ---\n")
cat("Probes remaining:", nrow(mSetSq), "\n")
table(keep)

# 10D: Remove cross-reactive probes (using maxprobes)

xreactive <- xreactive_probes(array_type = "450K")
keep <- !(featureNames(mSetSq) %in% xreactive)
mSetSq <- mSetSq[keep, ]
cat("\n--- After cross-reactive probe removal ---\n")
cat("Cross-reactive probes removed:", sum(!keep), "\n")
cat("Probes remaining:", nrow(mSetSq), "\n")
table(keep)

# 10E: Remove polymorphic CpGs (using minfi's dropLociWithSnps)
# This removes probes where SNPs are at the CpG interrogation site or 
# single base extension site, with minor allele frequency > 1%
mSetSq <- dropLociWithSnps(mSetSq, snps = c("CpG", "SBE"), maf = 0.01)
cat("\n--- After polymorphic CpG removal ---\n")
cat("Probes remaining:", nrow(mSetSq), "\n")

# Final summary
cat("\n", rep("=", 60), "\n", sep = "")
cat("FINAL CLEAN DATASET SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")
cat("Starting probes (450K array):        485,512\n")
cat("Final probes after all filters:     ", nrow(mSetSq), "\n")
cat("Total probes removed:               ", 485512 - nrow(mSetSq), "\n")
cat("Percentage retained:                ", round(100 * nrow(mSetSq) / 485512, 1), "%\n")
cat(rep("=", 60), "\n\n", sep = "")
print(mSetSq)

# ---- SECTION 11: VISUALIZE BEFORE vs AFTER NORMALIZATION --------------------
par(mfrow = c(1, 2))

densityPlot(rgset,
            sampGroups = pData(rgset)$Group,
            main = "Raw", legend = FALSE)
legend("top", 
       legend = levels(factor(pData(rgset)$Group)),
       text.col = brewer.pal(8, "Dark2"))

densityPlot(getBeta(mSetSq),
            sampGroups = pData(rgset)$Group,
            main = "Normalized", legend = FALSE)
legend("top", 
       legend = levels(factor(pData(rgset)$Group)),
       text.col = brewer.pal(8, "Dark2"))

dev.print(png, "Raw_vs_Normalized.png", width = 1200, height = 600, res = 150)

cat("Saved: Raw_vs_Normalized.png\n")

# ---- SECTION 12: EXTRACT BETA AND M-VALUES ----------------------------------
beta <- getBeta(mSetSq)
Mval <- getM(mSetSq)

cat("\n--- Beta value matrix dimensions ---\n")
print(dim(beta))
cat("\n--- M-value matrix dimensions ---\n")
print(dim(Mval))
cat("\n--- Head of Beta values (first 5 samples, first 5 probes) ---\n")
print(head(beta[, 1:5]))
cat("\n--- Head of M values (first 5 samples, first 5 probes) ---\n")
print(head(Mval[, 1:5]))

# ---- SECTION 13: Save Beta and M-value density plots ----
par(mfrow = c(1, 2))

densityPlot(beta,
            sampGroups = pData(mSetSq)$Group,
            main = "Beta values",
            legend = FALSE,
            xlab = "Beta values")
legend("top",
       legend = levels(factor(pData(mSetSq)$Group)),
       text.col = brewer.pal(8, "Dark2"),
       cex = 0.8)

densityPlot(Mval,
            sampGroups = pData(mSetSq)$Group,
            main = "M-values",
            legend = FALSE,
            xlab = "M values")
legend("topleft",
       legend = levels(factor(pData(mSetSq)$Group)),
       text.col = brewer.pal(8, "Dark2"),
       cex = 0.8)

dev.print(png, "BetaMvalue_Distributions.png", width = 1200, height = 600, res = 150)

cat("Saved: BetaMvalue_Distributions.png\n")

# ---- SECTION 14: MDS Plots ----
pal <- brewer.pal(8, "Dark2")
par(mfrow = c(1, 2))

plotMDS(getM(mSetSq),
        top = 1000,
        gene.selection = "common",
        col = pal[factor(pData(mSetSq)$Group)],
        labels = NULL,
        pch = 16,
        main = "MDS (dim 1 vs 2)")
legend("topright",
       legend = levels(factor(pData(mSetSq)$Group)),
       col = pal[1:2],
       pch = 16,
       cex = 0.8)

plotMDS(getM(mSetSq),
        top = 1000,
        gene.selection = "common",
        col = pal[factor(pData(mSetSq)$Group)],
        labels = NULL,
        pch = 16,
        dim = c(1, 3),
        main = "MDS (dim 1 vs 3)")
legend("topright",
       legend = levels(factor(pData(mSetSq)$Group)),
       col = pal[1:2],
       pch = 16,
       cex = 0.8)

dev.print(png, "MDS_Plots.png", width = 1400, height = 700, res = 150)

cat("Saved: MDS_Plots.png\n")
# ---- SECTION 15: DIFFERENTIAL METHYLATION — PROBE-LEVEL (DMPs) --------------

celltype <- factor(pData(mSetSq)$Group)
design   <- model.matrix(~0 + celltype)
colnames(design) <- levels(celltype)
design

# Fit linear model on M-values
fit <- lmFit(Mval, design) 

# Define contrast: DCIS vs adjacent-normal
contMatrix <- makeContrasts(
  DCIS_vs_Normal = DCIS - Adjacent.Normal,
  levels         = design)

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))

# ---- SECTION 16: ANNOTATE AND EXTRACT DMP RESULTS ---------------------------
ann450kSub <- anno450k[match(rownames(Mval), anno450k$Name),
                       c(1:4, 12:19, 24:ncol(anno450k))]
DMPs <- topTable(
  fit2,
  number   = Inf,
  coef     = "DCIS_vs_Normal",
  genelist = ann450kSub
)
head(DMPs)

cat("Total probes in results:", nrow(DMPs), "\n")

write.csv(DMPs, file = "DMPs_DCIS_vs_Normal.csv", row.names = FALSE)
cat("Saved: DMPs_DCIS_vs_Normal.csv\n")

# ---- SECTION 17: TOP CpG PLOTS (differentially methylated CpGs) -----------------------------------------------

par(mfrow = c(2, 2))
sapply(rownames(DMPs)[1:4], function(cpg) {
  plotCpg(beta, cpg = cpg, pheno = pData(mSetSq)$Group, ylab = "Beta values")
})

dev.print(png, "Top4_CpG_Plots.png", width = 1200, height = 1000, res = 150)
cat("Saved: Top4_CpG_Plots.png\n")

# ---- SECTION 18: VOLCANO PLOTS -----------------------------------------------

# 18a. Basic ggplot2 volcano

volcano_data <- data.frame(
  logFC = DMPs$logFC,
  negLogP = -log10(DMPs$P.Value),
  Probe = rownames(DMPs)
)

p_volcano <- ggplot(volcano_data, aes(x = logFC, y = negLogP)) +
  geom_point(alpha = 0.4, size = 0.5, color = "steelblue") +
  theme_minimal() +
  labs(
    title = "Volcano Plot - DCIS vs Adjacent-Normal",   
    x     = "Log Fold Change (logFC)",
    y     = "-log10(p-value)"
  ) +
  geom_hline(yintercept = -log10(0.05), col = "red",       linetype = "dashed") +
  geom_vline(xintercept = c(-0.2, 0.2),  col = "darkgreen", linetype = "dashed")

# Save the ggplot volcano
ggsave("VolcanoPlot_ggplot2.png", width = 6, height = 5, dpi = 300)

# 18b. EnhancedVolcano 

png("EnhancedVolcano_DMPs.png", width = 1400, height = 1200, res = 150)
EnhancedVolcano(
  DMPs,
  lab              = rownames(DMPs),
  x                = "logFC",      # column in DMps for fold change
  y                = "adj.P.Val",  # column in DMPs for significance
  pCutoff          = 0.05,
  FCcutoff         = 0.2,
  pointSize        = 1.5,
  labSize          = 2.5,
  title            = "DCIS vs Adjacent-Normal",
  subtitle         = "Differentially Methylated Probes (450K)",
  legendPosition   = "right",
  colAlpha         = 0.7,
  gridlines.major  = FALSE,
  gridlines.minor  = FALSE,
  caption          = paste0("Total probes = ", nrow(DMPs))
)
dev.off()
cat("Saved: EnhancedVolcano_DMPs.png\n")


# ---- SECTION 19: ANNOTATE DMPs WITH GENOMIC CONTEXT -------------------------
# Only adds annotation columns not already present in DMPs (avoids duplicates)

extra_cols <- c("chr", "pos", "UCSC_RefGene_Name",
                "UCSC_RefGene_Group", "Relation_to_Island")
new_cols   <- setdiff(extra_cols, colnames(DMPs))

if (length(new_cols) > 0) {
  DMPs_annot <- cbind(DMPs,
                      anno450k[match(rownames(DMPs), anno450k$Name),
                               new_cols, drop = FALSE])
} else {
  DMPs_annot <- DMPs
  cat("Note: All annotation columns already present - no duplication.\n")}

write.csv(DMPs_annot, file = "DMPs_DCIS_vs_Normal_annotated.csv", row.names = TRUE)
cat("Saved: DMPs_DCIS_vs_Normal_annotated.csv\n")

# Significant DMP overview
sigDMPs <- subset(DMPs_annot, adj.P.Val < 0.05 & abs(logFC) > 0.2)

cat("\n--- Significant DMP count ---\n")
print(nrow(sigDMPs))

cat("\n--- CpG island relation breakdown (significant DMPs) ---\n")
print(table(sigDMPs$Relation_to_Island))

cat("\n--- Top 10 gene region categories (significant DMPs) ---\n")
# Extracts only the first listed region per probe to avoid the hundreds of
# semicolon-concatenated combinations that table() produces otherwise
first_region <- sapply(strsplit(sigDMPs$UCSC_RefGene_Group, ";"), `[`, 1)
print(sort(table(first_region), decreasing = TRUE)[1:10])

# ---- SECTION 20: DMR ANALYSIS (Differentially Methylated Regions) -----------

myAnnotation <- cpg.annotate(
  object        = Mval,
  datatype      = "array",
  what          = "M",
  analysis.type = "differential",
  design        = design,
  contrasts     = TRUE,
  cont.matrix   = contMatrix,
  coef          = "DCIS_vs_Normal",   # must match contrast name
  arraytype     = "450K"
)


str(myAnnotation)


DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)

cat("\n--- Top DMRs ---\n")
print(results.ranges)

# Save DMR results
write.csv(as.data.frame(results.ranges),
          file = "DMRs_DCIS_vs_Normal.csv",
          row.names = FALSE)
cat("Saved: DMRs_DCIS_vs_Normal.csv\n")



# ---- SECTION 21: GENE ONTOLOGY TESTING (missMethyl) -------------------------

sigCpGs <- rownames(subset(DMPs, adj.P.Val < 0.05 & abs(logFC) > 0.2))
allCpGs <- rownames(DMPs)

cat("\n--- Number of significant CpGs for GO testing ---\n")
print(length(sigCpGs))

# gometh corrects for probe-number bias per gene
png("GO_BiasCorrectionPlot.png", width = 800, height = 600, res = 150)
gsa <- gometh(
  sig.cpg   = sigCpGs,
  all.cpg   = allCpGs,
  plot.bias = TRUE
)
dev.off()
#cat("Saved: GO_BiasCorrectionPlot.png\\n")
cat("\n--- Top 10 enriched GO terms ---\n")
print(topGSA(gsa, number = 10))

write.csv(topGSA(gsa, number = 50),
          file = "GO_enrichment_DCIS_vs_Normal.csv",
          row.names = FALSE)
cat("Saved: GO_enrichment_DCIS_vs_Normal.csv\n")



# ---- SECTION 22: SAVE SESSION INFO ------------------------------------------

tryCatch({
  si <- suppressWarnings(sessionInfo())
  writeLines(capture.output(si), "sessionInfo.txt")
  cat("Saved: sessionInfo.txt\n")
}, error = function(e) {
  cat("Warning: sessionInfo() failed. Saving fallback version info.\n")
  si_lines <- c(
    paste("R version:", R.version$version.string),
    paste("Date:", Sys.time()),
    "",
    "Loaded namespaces:",
    paste(sort(loadedNamespaces()), collapse = ", ")
  )
  writeLines(si_lines, "sessionInfo.txt")
  cat("Saved: sessionInfo.txt (fallback version)\n")
})


# How many DMRs were found?
cat("Number of DMRs found:", length(results.ranges), "\n")

# View the top 10 DMRs
print(head(as.data.frame(results.ranges), 10))

cat("\n=== BIOLOGICAL SUMMARY ===\n")
cat("Hypermethylated probes in DCIS:", sum(sigDMPs$logFC > 0), "\n")
cat("Hypomethylated probes in DCIS:", sum(sigDMPs$logFC < 0), "\n")



# ---- EXPORT ML FEATURE MATRIX (VARIANCE-PREFILTERED, LEAK-FREE) ------------

probe_var <- apply(beta, 1, var)
n_prefilter <- min(20000, length(probe_var))
top_var_cpgs <- names(sort(probe_var, decreasing = TRUE))[1:n_prefilter]

beta_subset <- beta[top_var_cpgs, ]
beta_subset_t <- t(beta_subset)

labels <- pData(mSetSq)$Group
beta_df <- as.data.frame(beta_subset_t)
beta_df$label <- ifelse(labels == "DCIS", 1, 0)

if (!dir.exists("results")) dir.create("results")
write.csv(beta_df,
          "results/beta_variance_prefiltered.csv",
          row.names = TRUE)
cat("Exported:", nrow(beta_df), "samples x", ncol(beta_df) - 1,
    "variance-prefiltered CpGs (label-blind selection)\n")
