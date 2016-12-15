clusterization = function(mat.exp,wayOut,name.data,mtClust,matDist){
  
  # mat.exp = GeNNet
  # wayOut = dirW
  # name.data = nameExp
  # mtClust = "median"
  # matDist = "Pearson"
  
  # --- Colors pallet ---
  # === colores to heatMap ===
  cols = colorRampPalette(brewer.pal(10, "RdBu"))(256)
  
  # =====================
  # --- Plot Heat Map ---
  # =====================
  expMat = exprs(mat.exp[[1]])
  expMat = expMat[which(is.na(match(rownames(expMat),mat.exp$selected$probes)) == FALSE),]
  
  rownames(expMat) = mat.exp$selected[which(is.na(match(rownames(expMat),mat.exp$selected$probes)) == FALSE),"symbol"]
  
  phenos = pData(mat.exp[[1]])  
  colNamesD = colnames(expMat)
  sets = phenos[which(is.na(match(colNamesD,phenos$SAMPLE_NAME))==F),"SETS"]
  colnames(expMat) = gsub(".CEL.gz","",colnames(expMat))
  
  rainbowcols = rainbow(length(unique(phenos$SETS)),s = 0.5)
  matColor = as.array(rainbowcols)
  names(matColor) = unique(phenos$SETS)
  
  # --- define colors in sets ---
  cow = NULL
  for(i in 1:length(sets)){
    cow =  c(cow,as.character(matColor[match(sets[i],names(matColor))] ))
  }
  patientcolors  = as.matrix(cow)
  colnames(patientcolors) = "sets"
  
  # --- Calculating z-scores for the rows ---
  data = expMat 
  data <- t(scale(t(data)))
  all.na <- apply(data,1, function(x) all(is.na(x)))
  data <- data[!all.na,]
  expMat = data
  
  # --- Metricas ---
  if(matDist == "Euclidian"){
    distance = dist(expMat,method="euclidian") 
  }
  if(matDist == "Manhattan"){
    distance = dist(expMat,method="manhattan") 
  }
  if(matDist == "Pearson"){
    dissimilarity = 1 - cor(t(expMat))
    distance = as.dist(dissimilarity)
  }
  tree.corr = hclust(distance, method=mtClust)
  
  # =====================
  # === Plot Heat Map ===
  # =====================
  source("Module-A/HeatMap3.R")
  pdf(paste(wayOut,name.data,"-HeatMap.pdf",sep="")) 
  main_title=paste(name.data,"- experiment")
  heatmap.3(expMat, 
            # hclustfun=myclust,
            Rowv=as.dendrogram(tree.corr),
            # distfun=mydist, 
            # na.rm = TRUE, 
            scale="none",#"row",
            dendrogram="both", 
            margins=c(6,12),
            Colv=TRUE, 
            ColSideColors=patientcolors,
            symbreaks=FALSE, 
            key=TRUE, 
            symkey=FALSE,
            density.info="none", 
            trace="none", 
            main=main_title, 
            # labCol=FALSE, 
            labRow="",
            cexRow=1,
            #density="none", 
            col=rev(cols),
            ColSideColorsSize=2, 
            RowSideColorsSize=2, 
            KeyValueName="z-score")
  legend("topright",legend=c(names(matColor)),
         fill=matColor, border=FALSE, bty="n", y.intersp = 0.9, cex=1)
  dev.off()
  
  # ===========
  # --- fpc ---
  # ===========
  options(digits=3)
  set.seed(20000)
  # face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
  pk1 = pamk(distance,krange=1:15,criterion="ch",critout=TRUE)
  pk1$nc
  # names(pk1)
  pk1$pamobject
  pk1$crit
  
  # --- BIND COLS fpc ---	
  groups= unlist(pk1[[1]]["clustering"])
  mat.fpc = as.data.frame(cbind(expMat,groups))
  
  name=paste(name.data,"-PAM",sep="")
  
  # plot_cluster(mat.fpc,wayOut,name)
  
  mat.fpc2 = as.array(mat.fpc[,"groups"])
  rownames(mat.fpc2) = rownames(mat.fpc)
  mat.fpc2
  
  names.list.Exp = c(names(mat.exp),"clusterization")
  mat.exp[[6]] = mat.fpc2
  names(mat.exp) = names.list.Exp  
  
  
  ClustPar = list(matExp=expMat,colores=patientcolors,tree=tree.corr,matColor=matColor)
  save(ClustPar,file = paste(wayOut,"ClutPar.rda",sep=""))
  
  mat.exp
}