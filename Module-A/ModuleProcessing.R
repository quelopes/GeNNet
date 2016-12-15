moduleProcessing = function(){
  
  # === PROVENANCE ===
  library(RDataTracker)
  ddg.init(enable.console=TRUE)

  # === START ===
  ddg.start("dependences")
  source("Module-A/Dependences.R")
  dependences()
  graph = startGraph("http://localhost:7474/db/data/",username="neo4j",password = "graph")
  ddg.finish("dependences")
  
  # === CREATE NODES EXPERIMENT ===
  ddg.start("node experiment")
  cNodeExperiment(infoExp,graph)
  ddg.finish("node experiment")
  
  # === NORMALIZED ===
  ddg.start("normalize")
  GeNNet = pre_Raw(dirR_raw,dirW,pheno,method,nameExp)
  ddg.finish("normalize")
  
  # === E-SET ===
  ddg.start("e-set")
  GeNNet = pre_eSet(GeNNet,pheno,orgAnnot,nameExp)
  ddg.finish("e-set")
  
  # === LIMMA ===
  ddg.start("limma")
  GeNNet = pre_Limma(GeNNet[[1]],FDR,logFC,dirW,names(GeNNet[1]))
  ddg.finish("limma")
  
  # === ANNOTATION ===
  ddg.start("annotation")
  GeNNet = pre_PosLimma(mat_annot,GeNNet,nameExp,dirW)
  ddg.finish("annotation")
  
  # === NODES GENES EXPRESSION ===
  ddg.start("nodes expression")
  cNodeValues(GeNNet,graph)
  ddg.finish("nodes expression")
  
  if(nrow(GeNNet$selected) >0){
    
    # === CLUSTERIZATION ===
    ddg.start("clusterization")
    GeNNet = clusterization(GeNNet,dirW,nameExp,mClust,matDist)
    ddg.finish("clusterization")
    
    # === NODES CLUSTERS ===
    ddg.start("nodes cluster")
    clust_pop(GeNNet,paste(mClust,matDist,sep="-"),graph)
    ddg.finish("nodes cluster")
    
    # === WGCNA ===
    ddg.start("nodes wgcna")
    power = powerDetec(GeNNet,nameExp,dirW)
    amd = amdDetec(power,GeNNet,nameExp,dirW)
    ddg.finish("nodes wgcna")
    
    # === ENRICHMENT ===  
    ddg.start("functional analysis")
    GeNNet = enrichment_GOstats2(GeNNet,nameExp,orgAnnot,dirW)
    ddg.finish("functional analysis")
    
    # === NODES ENRICHMENT ===
    ddg.start("nodes funct.analysis")
    pop_enrichment(GeNNet,paste(mClust,matDist,sep="-"),graph)
    ddg.finish("nodes funct.analysis")
  }
  
  # === DATA PERSISTENCE AND STATISTICS ===
  ddg.start("persistence")
  GeNNet = graphDataPersist(GeNNet,graph,dirW)  
  
  # === SAVE RDATA ===  
  save(GeNNet,file=paste(dirW,nameExp,"-GeNNet.rda",sep=""))
  ddg.finish("persistence")
  
  ddg.save()
  file.copy(from="/home/rstudio/ddg/", to=dirW, 
            overwrite = TRUE, recursive = TRUE, 
            copy.mode = TRUE)
}