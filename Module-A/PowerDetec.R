powerDetec = function(mat,n.exp,out){
  
  # mat = GeNNet
  # n.exp = nameExp  
  # out = wayOut
  
  
  expMat = exprs(mat[[1]])
  expMat = expMat[which(is.na(match(rownames(expMat),mat$selected$probes)) == FALSE),]
  
  rownames(expMat) = mat$selected[which(is.na(match(rownames(expMat),mat$selected$probes)) == FALSE),"symbol"]
  
  # --- define sequences ---
  powers = c(1:15,seq(16,30,2))
  
  sft = pickSoftThreshold(t(expMat), powerVector = powers, verbose = 5)
  
  powerRanges = sft$powerEstimate
  Power = sft$fitIndices
  
  # ====================
  # === Plot figures ===
  # ====================
  pdf(paste(out,"Soft-Threshold-",n.exp,".pdf",sep=""))
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence - ",n.exp,sep=""));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity - ",n.exp,sep=""))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  # par(mfrow=c(1,1))
  
  
  pst = sft$fitIndices
  # ====================
  # === Plot figures ===
  # ====================
  pdf(paste(out,"Soft-Threshold-Values-",n.exp,".pdf",sep=""))
  par(mfrow=c(2,2), mar=c(4,4,3,2)+0.1)
  plot(powers, -sign(pst[,3])*pst[,2],xlab="Soft Threshold (power)", 
       ylab="SFT Model Fit, signed R^2", type="n", main=paste("Scale Independence - ",n.exp,sep=""))
  text(pst[,1], -sign(pst[,3])*pst[,2],labels=powers,col="red")
  abline(h=0.37,col="red")
  
  plot(pst[,1], pst$Density, type="n",xlab="Soft Threshold (power)",
       ylab="Density", main=paste("Density - ",n.exp,sep=""))
  text(pst[,1],pst$Density,labels=powers,col="red")
  
  plot(pst[,1], pst$Heterogeneity, type="n",xlab="Soft Threshold (power)",
       ylab="Heterogeneity", main=paste("Heterogeneity - ",n.exp,sep=""))
  text(pst[,1],pst$Heterogeneity,labels=powers,col="red")
  
  plot(pst[,1], pst$Centralization, type="n",xlab="Soft Threshold (power)",
       ylab="Centralization", main=paste("Centralization - ",n.exp,sep=""))
  text(pst[,1],pst$Centralization,labels=powers,col="red")
  
  dev.off()
  # par(mfrow=c(1,1))
  
  
  powerRanges
  
}