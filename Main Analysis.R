# ==============================================================================
# Project: Machine Learning-Driven Discovery of VEGFC as a Prognostic Cytokine 
# Biomarker in Breast Cancer
# Script: main_analysis.R
# Description: This script contains the core computational pipeline, including 
# LASSO feature selection, survival analysis, immune landscape characterization, 
# and pharmacogenomic screening.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Packages
# ------------------------------------------------------------------------------
library(glmnet)      # For LASSO machine learning
library(survival)    # For survival analysis
library(survminer)   # For Kaplan-Meier plots
library(ggplot2)     # For data visualization
library(pROC)        # For ROC curve analysis
library(TCGAplot)    # For pan-cancer and clinical association plots
library(oncoPredict) # For pharmacogenomic screening
library(clusterProfiler) # For GSEA analysis

# ------------------------------------------------------------------------------
# 1. Data Acquisition and Processing 
# ------------------------------------------------------------------------------
# Load TCGA-BRCA RNA-seq and clinical data
# Note: Replace 'path/to/your/data.csv' with your actual local file paths
tcga_expr <- read.csv("data/TCGA_BRCA_expression.csv", row.names = 1)
tcga_clin <- read.csv("data/TCGA_BRCA_clinical.csv", row.names = 1)

# Filter for relevant cytokine panel genes
cytokine_panel <- c("VEGFC", "IL6", "TNF", "CCL2", "CXCL12") # Add full panel here
tcga_expr_filtered <- tcga_expr[rownames(tcga_expr) %in% cytokine_panel, ]

# ------------------------------------------------------------------------------
# 2. Machine Learning-Driven Feature Selection (LASSO) 
# ------------------------------------------------------------------------------
# Prepare matrix for glmnet
x_matrix <- t(tcga_expr_filtered)
# Define response variable: Nodal status (N0 = 0, non-N0 = 1)
y_response <- as.numeric(tcga_clin$N_stage != "N0") 

# Perform 10-fold cross-validation to find optimal lambda
set.seed(123)
cv_fit <- cv.glmnet(x_matrix, y_response, family = "binomial", alpha = 1)

# Plot Cross-Validation curve 
plot(cv_fit)

# Extract coefficients at optimal lambda (lambda.min)
lasso_coef <- coef(cv_fit, s = "lambda.min")
print(lasso_coef) # VEGFC should emerge as the most significant feature

# ------------------------------------------------------------------------------
# 3. Clinical Association and Survival Analysis (Section 2.3)
# ------------------------------------------------------------------------------
# ROC Curve Analysis 
roc_obj <- roc(tcga_clin$Tissue_Type, tcga_expr["VEGFC", ])
plot(roc_obj, print.auc = TRUE)

# Relapse-free survival (RFS) analysis 
# Stratify by median VEGFC expression
median_vegfc <- median(as.numeric(tcga_expr["VEGFC", ]))
tcga_clin$VEGFC_Group <- ifelse(tcga_expr["VEGFC", ] > median_vegfc, "High", "Low")

# Kaplan-Meier Curve
km_fit <- survfit(Surv(RFS_time, RFS_status) ~ VEGFC_Group, data = tcga_clin)
ggsurvplot(km_fit, data = tcga_clin, pval = TRUE, risk.table = TRUE)

# Multivariate Cox proportional hazards model 
cox_model <- coxph(Surv(RFS_time, RFS_status) ~ VEGFC_Group + Age + Stage + ER_status + HER2_status, data = tcga_clin)
summary(cox_model)
ggforest(cox_model, data = tcga_clin)

# ------------------------------------------------------------------------------
# 4. Immune Microenvironment and Pan-Cancer Analysis 
# ------------------------------------------------------------------------------
# Note: TCGAplot functions are used here as described in the methodology
# Example: Pan-cancer correlation of VEGFC with immune checkpoints 
# TCGAplot::pan_cancer_cor(gene1 = "VEGFC", gene_list = c("PDCD1", "CTLA4", "TIGIT", "HAVCR2"))

# ------------------------------------------------------------------------------
# 5. Pharmacogenomic Screening 
# ------------------------------------------------------------------------------
# Predict drug sensitivity (IC50) using oncoPredict
# calcPhenotype(trainingExprData = ..., trainingPtype = ..., testExprData = as.matrix(tcga_expr))
# Correlate predicted IC50 of AZD6482 with VEGFC expression
# cor.test(tcga_expr["VEGFC", ], predicted_IC50_AZD6482, method = "pearson")

# End of Script