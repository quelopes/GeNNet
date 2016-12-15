amdDetec = function(pwCalc,mat,n.exp,out){
  
  # out = wayOut
  # mat = GeNNet
  # pwCalc = power
  # n.exp = nameExp
  
  
  expMat = exprs(mat[[1]])
  expMat = expMat[which(is.na(match(rownames(expMat),mat$selected$probes)) == FALSE),]
  
  rownames(expMat) = mat$selected[which(is.na(match(rownames(expMat),mat$selected$probes)) == FALSE),"symbol"]
  
  # ------------------
  # --- Parameters ---
  # ------------------
  minModuleSize = 30
  maxBlockSize=1500
  corType="pearson"
  mergingThresh = 0.25
  
  mat = t(expMat)
  amd = list()
  if(is.na(pwCalc)==TRUE){
    pwCalc = 9
  }
  
  net = blockwiseModules(mat,corType=corType,
                         maxBlockSize=maxBlockSize,networkType="signed",power=pwCalc,minModuleSize=minModuleSize,
                         mergeCutHeight=mergingThresh,numericLabels=TRUE,saveTOMs=TRUE,TOMType = "signed",
                         pamRespectsDendro=FALSE,saveTOMFileBase=paste(out,"grupo-",n.exp,"-TOM",sep=""))
  mergedColors = labels2colors(net$colors)
  amd = list(data = mat, colors = mergedColors, net = net)
  
  # ver exemplos para usar depois
  moduleLabels = net$colors
  #moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  
  
  # -------------------
  # --- Plot figure ---
  # -------------------
  # pdf(paste(out,"Dendrogram-",n.exp,".pdf",sep=""))
  # plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],main = paste("Cluster Dendrogram ",n.exp,sep=""),
  #                     "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
  # dev.off()
  
  # -------------------
  # --- Plot figure ---
  # -------------------
  table(net$colors)
  #sizeGrWindow(12,9)
  # mergedColors = labels2colors(net$colors)
  pdf(paste(out,"Plot_dendro_and_colors-AMD-",n.exp,".pdf",sep=""))
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels= FALSE, hang=0.03, addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  
  # --------------
  # --- saving ---
  # --------------
  
  
  names(amd) = c(n.exp,names(amd[2]),names(amd[3]))
  
  #   geneTree = net$dendrograms[[1]]
  
  save(amd, file = paste(out,"AMD_results.RData",sep=""))
  
  names(net)
  # [1] "colors"         "unmergedColors" "MEs"            "goodSamples"    "goodGenes"      "dendrograms"    "TOMFiles"      
  # [8] "blockGenes"     "blocks"         "MEsOK
  
  amd
  
}