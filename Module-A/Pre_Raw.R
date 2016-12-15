pre_Raw = function(wayRaw,wayOut,pheno,method,nameExp){
  
  # wayRaw = dirR_raw
  # wayOut = dirW
  # method = method
  # nameExp = nameExp
  # pheno
  
  print("*** ACTIVITY NORMALIZATION ***")
  
  goAfter = getwd()
  setwd(wayRaw)
  raw = ReadAffy()
  
  # === PLOTS ===
  
  # --- Pallet ---
  groupnames = data.frame(sets=pheno$SETS)
  # tmData = unique(groupnames$sets)
  tmData = names(table(groupnames$sets))
  if(length(tmData) == 2){
    cols = brewer.pal(3, "Dark2")
  }
  if(length(tmData) > 2){
    cols = brewer.pal(tmData, "Dark2")
  }
  groupnames = cbind(groupnames,colores=NA)
  for(i in 1:length(tmData)){
    sel = which(groupnames$sets==tmData[i])
    groupnames[sel,"colores"] = cols[i]
  }
  
  # --- Boxplot of probe level data for each array ---
  pdf(paste(wayOut,nameExp,"Quality-Before-Norm-Box-plots-Of-Probe-Level-Data",".pdf",sep=""),family="Helvetica", width=640, height=480, paper="a4r")
  boxplot(raw, main="Probe level data", col=groupnames$colores, cex.axis=0.75)
  axis(side=1, line=2, 
       # at=c(nsBLC/2, nsBLC+(nsnBLC/2)), 
       # labels=tmData, 
       tick=FALSE, 
       font=2, 
       cex=1.2)
  legend("topright",tmData,fil = unique(groupnames$colores))
  dev.off()
  
  # --- Density estimates ---
  pdf(paste(wayOut,nameExp,"Quality-before-norm-Density-estimates-of-probe-level-data",".pdf",sep=""),family="Helvetica", width=640, height=480, paper="a4r")
  hist(raw, col=cols, xlab="log2 intensities")
  dev.off()
  
  # --- Spearman correlation ---
  raw.data2 = exprs(raw)      
  # Compute and view concordance correlation and Spearman correlation (especially the computation of the concordance correlation takes time)
  spearman.cor = cor(raw.data2,  method=c("spearman"))
  # round(spearman.cor, 3)
  # mean(spearman.cor)                            # mean correlation over all arrays
  # round(apply(spearman.cor, 2, median), 3)      # median correlation per array
  
  pdf(paste(wayOut,nameExp,"Quality-before-norm-Spearman-correlation-of-probe-level-data",".pdf",sep=""),family="Helvetica",paper="a4r")
  heatmap(x=-spearman.cor, main="Spearman correlation of probe level data", scale="none", col=heat.colors(256))
  dev.off()
  
  
  # === NORMALIZATION ===
  # --- if mas5 --- 
  if(method == "mas5"){
    eset.mas5 = mas5(raw) # few minutes...
    exprSet.nologs = exprs(eset.mas5)
    exprSet = log(exprSet.nologs, 2)
    #summary(exprSet)
  }
  
  # --- if rma --- 
  if(method == "rma"){
    eset.rma = rma(raw) # few minutes...
    exprSet.nologs = exprs(eset.rma)
    exprSet = log(exprSet.nologs, 2)
    #summary(exprSet)
  }
  
  # --- if gcrma --- 
  if(method == "gcrma"){
    raw = ReadAffy()
    eset.gcrma = gcrma(raw)
  }
  
  # --- if mas5calls --- 
  if(method == "mas5calls"){
    raw = ReadAffy()
    eset.mas5calls = mas5calls(raw)
  }
  
  # --- if justPlier --- 
  if(method == "justPlier"){
    raw = ReadAffy()
    eset.justPlier = justPlier(raw)
  }
  
  # --------------------------------
  # --- Plot after normalization ---
  # --------------------------------
  pdf(paste(wayOut,nameExp,"Quality-After-Norm-Box-plots-Of-Probe-Level-Data",".pdf",sep=""),family="Helvetica", width=640, height=480, paper="a4r")
  boxplot(exprSet, main="Probe level data", col=groupnames$colores, cex.axis=0.75)
  axis(side=1, line=2, 
       # at=c(nsBLC/2, nsBLC+(nsnBLC/2)), 
       # labels=tmData, 
       tick=FALSE, 
       font=2, 
       cex=1.2)
  legend("topright",tmData,fil = unique(groupnames$colores))
  dev.off()
  
  # --- Density estimates ---
  pdf(paste(wayOut,nameExp,"Quality-After-norm-Density-estimates-of-probe-level-data",".pdf",sep=""),family="Helvetica", width=640, height=480, paper="a4r")
  hist(exprSet, xlab="log2 intensities")#col=cols, 
  dev.off()
  
  # --- Spearman correlation ---
  raw.data2 = exprSet      
  # Compute and view concordance correlation and Spearman correlation (especially the computation of the concordance correlation takes time)
  spearman.cor = cor(raw.data2,  method=c("spearman"))
  # round(spearman.cor, 3)
  # mean(spearman.cor)                            # mean correlation over all arrays
  # round(apply(spearman.cor, 2, median), 3)      # median correlation per array
  
  pdf(paste(wayOut,nameExp,"Quality-After-norm-Spearman-correlation-of-probe-level-data",".pdf",sep=""),family="Helvetica",paper="a4r")
  heatmap(x=-spearman.cor, main="Spearman correlation of probe level data", scale="none", col=heat.colors(256),cexRow=0.7,cexCol=0.7)
  dev.off()
  
  setwd(goAfter)
  
  # --- writing data (table) ---
  write.table(exprSet, file=paste(wayOut,nameExp,"-MatrixNormalized.csv",sep=""), quote=F, sep=",")
  exprSet
}
