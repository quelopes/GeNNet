dependences = function(){
  
  # ================================
  # ---------- PACKAGES ------------
  # ================================
  
  print("*** STARTING LIBRARIES AND FUNCTIONS ***")	
  
  # ----- PRE PROCESSING -----	
  library(affy) # Pre_Raw
  library(impute) # Pre_Raw
  library(Biobase) # S4 Read Esets
  library(limma) # Pre_Limma
  library(hgu133plus2cdf)
  library(hgu133a2cdf)
  library(hugene10stv1cdf)
  library(frma)
  
  # ----- ANNOTATION & ENRICHMENT -----
  library(org.Hs.eg.db) 
  library(org.Mmu.eg.db) 
  library(org.Mm.eg.db) 
  library(org.Rn.eg.db)
  
  # ----- PLOTS -----
  library(ggplot2) 
  library(gplots)
  library(RColorBrewer)
  
  # ----- GRAPH DATABASE -----
  library(RNeo4j)
  library(RDataTracker)
  library(igraph)
  
  # --- ANALYSIS ---
  library(igraph)
  library(VennDiagram)
  library(gplots) # Heatmap.2
  library(fpc) # cluster PAM
  library(stringr) # str_split
  library(GOstats) # Enrichment
  library(WGCNA)
  library(dynamicTreeCut)
  
  # ================================
  # --------- FUNCTIONS ------------
  # ================================
  
  # ----- PRE PROCESSING -----
  source("Module-A/Pre_Raw.R") 
  source("Module-A/Pre_eSet.R") 
  source("Module-A/Pre_Limma.R") 
  source("Module-A/Pre_PosLimma.R") 
  
  source("Module-A/MyStandardise.R")
  source("Module-A/Clusterization.R")
  source("Module-A/HeatMap3.R")
  source("Module-A/Enrichment_GOStat2.R")
  source("Module-A/PowerDetec.R")
  source("Module-A/AmdDetec.R")
  
  # ----- LOAD DATABASE & GET INFORMATION -----
  source("Module-A/cNodeExperiment.R")
  source("Module-A/cNodeValues.R")
  source("Module-A/Clust_pop.R") # pop Graph Database with nodes cluster
  source("Module-A/Pop_enrichment.R")  
  source("Module-A/GraphDataPersist.R")  
  
}