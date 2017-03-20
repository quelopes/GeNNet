# ======================
# === Configurations ===
# ======================
nameExp = "GSE62232"
orgAnnot = "org.Hs.eg.db"
GPL = "GPL570.txt"

# --- Experiment Information ---
overalDesign = "HCC liver tumors corresponding to 81 patients. In all cases, tumor samples were frozen (-80°C) after hepatic resection at diagnosis. Normal liver samples were used as reference samples. Comparative gene expression analysis was done using Affymetrix U133plus v2 array (GPL570)"
platName = "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array"
platAccess = "GPL570"

infoExp = list(overalDesign=overalDesign,platName=platName,platAccess=platAccess,namesExp=nameExp)

# --- write the results ---
dirW = paste(getwd(),"/Results/",sep="")

# --- normalization ---
method = "mas5"
# or "rma"
dirR_raw = paste("Data/CEL/",sep="")

# --- e-set ---
pheno = read.table(paste("Data/",nameExp,".csv",sep=""),head=T,sep=",")

# --- annotation ---
mat_annot = read.table(paste("Annotation/",GPL,sep=""),header=T,sep="\t",quote = "",fill=TRUE)

# ======================
# === GERAL SETTING === 
# ======================
# --- limma ---
FDR = 0.05 #threshold 
logFC = 1  #threshold

# --- clusterization ---
mClust = "average"
#"complete","average","median","centroid"

matDist = "Pearson"
# Spearman # Euclidian

source("Module-A/ModuleProcessing.R")
moduleProcessing()

