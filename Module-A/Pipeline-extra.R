#############################################################################################################
###                                                                                                       ###
### ACCOMPANYING COMPLETE R CODE FOR                                                                      ###
###                                                                                                       ###
### 'BIOINFORMATICS FOR OMICS DATA: METHODS AND PROTOCOLS'                                                ###
### CHAPTER 1.6: STATISTICS FOR OMICS DATA                                                                ###
### BY DANIELA DUNKLER (1), FÁTIMA SÁNCHEZ-CABO (2) AND GEORG HEINZE (1)                                  ###
###                                                                                                       ###
### (1) Medical University of Vienna; Center for Medical Statistics, Informatics and Intelligent Systems; ###
###     Section for Clinical Biometrics; Vienna; Austria                                                  ###
### (2) Centro Nacional de Investigaciones Cardiovasculares; Genomics Unit; Madrid; Spain                 ###
###                                                                                                       ###
### March 2010                                                                                            ###
###                                                                                                       ###
#############################################################################################################


#-----------------------------------------------------------------------------------------------------------#
#   Load required R packages                                                                                #
#-----------------------------------------------------------------------------------------------------------#
library(affy)                                              # read CEL files into R & preprocessing
library(affyPLM)                                           # fit robust linear model to probe level data
library(genefilter)                                        # non-specific filtering
library(limma)                                             # gene expression analysis
library(samr)                                              # SAM analysis
library(epiR)                                              # concordance correlation
library(mpm)                                               # spectral map
library(RankProd)                                          # RankProduct statistic
library(simpleaffy)                                        # Affymetrix quality assessment
library(RColorBrewer)                                      # colors for plots
library(grDevices)                                         # colors for plots
library(gdata)                                             # ‘keep’ function

# source("https://bioconductor.org/biocLite.R")
# biocLite(c("affyPLM","samr","epiR","mpm","simpleaffy","RankProd"))

#----------------------------------------------------------------------------------------------------------#
#    Set parameters manually                                                                               #
#----------------------------------------------------------------------------------------------------------#
phenotype <- c(rep("BLC", 8),rep("nBLC", 8))               # vector with the group membership of each array 
# (Arrays are sorted according to their name.)
input.folder  <- "/home/quelopes/DATA-home/GennetD/GSE3744/Data/CEL/"                                 # input folder with CEL files
output.folder <- "/home/quelopes/DATA-home/GennetD/GSE3744/Results/"                       # output folder where to save results

#----------------------------------------------------------------------------------------------------------#
#    Some functions                                                                                        #
#----------------------------------------------------------------------------------------------------------#

plot.cor.boxplot <- function(plot.data, plot.title, plot.cols, n.g1, n.g2, plot.groupnames)
{
  ####################################################################################
  ###                                                                              ###
  ### FUNCTION: Short function for plotting Boxplots of a correlation matrix       ###
  ###                                                                              ###
  ### INPUT:                                                                       ###
  ###   plot.data       = correlation matrix                                       ###
  ###   plot.title      = title of the plot                                        ###
  ###   plot.cols       = names of 2 colours                                       ###
  ###   n.g1            = number of samples in group 1                             ###
  ###   n.g2            = number of samples in group 2                             ###
  ###   plot.groupnames = groupnames                                               ###
  ###                                                                              ###
  ####################################################################################
  
  boxplot(data.frame(plot.data), main=plot.title, cex.axis=0.75, col=c(rep(plot.cols[1], n.g1), rep(plot.cols[2], n.g2)))
  abline(h=median(plot.data, na.rm=TRUE), col="gray", lty=2, lwd=2)
  mtext(text="Reference line is the median correlation", side=3, line=0.5, col="gray", font=2)
  axis(side=1, line=2, at=c(n.g1/2, n.g1+(n.g2/2)), tick=FALSE, font=2, cex=1.2, labels=c(plot.groupnames[1], plot.groupnames[2]))
}

genelist.sam <- function(x, y, fdr=0.01, nperms=100, t.sam=NA, delta.table=NA, refine.delta=TRUE, ...)
{
  ####################################################################################
  ###                                                                              ###
  ### FUNCTION: compute SAM analysis and save results                              ###
  ###                                                                              ###
  ### genelist.sam can be re-used with existing t.sam or delta.table objects       ###
  ###                                                                              ###
  ### INPUT:                                                                       ###
  ###   x            = filtered gene expression data (exprs(filtered.data))        ###
  ###   y            = group membership (c(rep(1,nsBLC),rep(2,nsnBLC)))            ###
  ###   fdr          = favored FDR (default 0.01)                                  ###
  ###   nperms       = number of permutations for SAM (default 100)                ###
  ###   t.sam        = object computed by samr function (optional)                 ###
  ###   delta.table  = delta table (optional)                                      ###
  ###   refine.delta = refine values for delta, TRUE or FALSE                      ###
  ###                                                                              ###
  ### OUTPUT: a list with genelist, t.sam, delta.table and pi0.                    ###
  ###                                                                              ###
  ####################################################################################
  
  data.sam <- list(x=x, geneid=as.character(1:nrow(x)),genenames=rownames(x), y=y,logged2=T)
  if(is.na(t.sam) | is.na(delta.table))
  {
    t.sam <- samr(data.sam, resp.type="Two class unpaired")
    if(!refine.delta) nvals=50
    else nvals=30
    delta.table <- samr.compute.delta.table(t.sam, nvals=nvals)      # use rough delta table to begin
  }
  
  delta.max <- delta.table[delta.table[,6] < fdr,1][1]               # select max delta for which FDR is controlled (90th percentile of FDR < FDRthreshold)
  delta.min <- delta.table[delta.table[,5] > 2*fdr, 1]               # select min delta for which FDR is controlled (median FDR > 2* FDRthreshold)
  delta.min <- as.vector(na.omit(delta.min))
  delta.min <- delta.min[length(delta.min)]
  # compute delta table once again, with more refined values for delta
  if(refine.delta == FALSE)
  {
    DEG.SAM <- samr.compute.siggenes.table(t.sam, delta.max, data.sam, delta.table, compute.localfdr=TRUE)
  }
  else
  {
    length.delta.table <- sum(-1*(is.na(delta.table[, 5]))+1)
    dels <- delta.min+(delta.max-delta.min)/length.delta.table*(0:length.delta.table)
    
    delta.table.fine <- samr.compute.delta.table(t.sam, dels=dels)    # use rough delta table to begin
    delta <- delta.table.fine[delta.table.fine[,5] < fdr,1][1]        # select delta for which FDR is controlled (median percentile of FDR < 0.01)
    DEG.SAM <- samr.compute.siggenes.table(t.sam, delta, data.sam, delta.table.fine, compute.localfdr=TRUE)
  }
  
  genelist <- rbind(DEG.SAM$genes.up, DEG.SAM$genes.lo)
  genelist[rank(-abs(as.numeric(genelist[,4]))),] <- genelist
  ngenes <- nrow(genelist)
  cat("Number of genes entered in analysis: ", nrow(x), "\n")
  cat("Estimated proportion of non-differentially expressed genes: ", t.sam$pi0, ".\n")
  cat("SAM found ", ngenes, " differentially expressed genes at an FDR of ", fdr, ".\n")
  list(genelist=genelist, t.sam=t.sam, delta.table=delta.table, pi0=t.sam$pi0, delta=delta, data.sam=data.sam)
}


plot.genes <- function(eset, group, geneid, qtable=NA, print.fc=FALSE, print.average=FALSE, print.sd=FALSE,
                       print.cindex=FALSE, print.qvalue=FALSE)
{
  ####################################################################################
  ###                                                                              ###
  ### FUNCTION: Plot gene expressions of a specified gene and optionally print     ###
  ###           some statistics                                                    ###
  ###                                                                              ###
  ### INPUT:                                                                       ###
  ###   eset          = matrix with gene expressions                               ###
  ###   group         = group membership                                           ###
  ###   genid         = which gene should be plotted                               ###
  ###   qtable        = q-table (required if you want to print some statistics)    ###
  ###   print.fc      = print fold change (TRUE or FALSE)                          ###
  ###   print.average = print average (TRUE or FALSE)                              ###
  ###   print.sd      = print standard deviation (TRUE or FALSE)                   ###
  ###   print.cindex  = print c-index (TRUE or FALSE)                              ###
  ###   print.qvalue  = print q-value (TRUE or FALSE)                              ###
  ###                                                                              ###
  ####################################################################################
  
  groupnames <- dimnames(table(group))
  groupnames <- c( groupnames[[1]][1], groupnames[[1]][2])
  n.group <- as.numeric(table(group))
  data.tmp <- eset[geneid,]
  
  plot(c(0,4), c(floor(min(data.tmp)), ceiling(max(data.tmp))), col="white", bty="n", xlab="", ylab="", main=geneid, axes=FALSE, cex=2)
  axis(1, at=c(1,3), labels=groupnames, cex=1)
  axis(2, cex=1)
  box()
  points(rep(1, n.group[1]), data.tmp[group==groupnames[1] ], pch=19)
  points(rep(3, n.group[2]), data.tmp[group==groupnames[2] ], pch=19)
  lines(c(0.7, 1.3), c(median(data.tmp[group==groupnames[1] ]), median(data.tmp[group==groupnames[1] ])))
  lines(c(2.7, 3.3), c(median(data.tmp[group==groupnames[2] ]), median(data.tmp[group==groupnames[2] ])))
  
  if(print.fc | print.average | print.sd | print.qvalue | print.cindex)
  {
    qvalue <- qtable[qtable[,1]==geneid,  6]
    average <- mean(data.tmp)
    sdev <- (sd(data.tmp[group==groupnames[1]])^2*n.group[1] + sd(data.tmp[group==groupnames[2]])^2*n.group[2])/sum(n.group)
    fc <- round(2^(mean(data.tmp[group==groupnames[1]]-mean(data.tmp[group==groupnames[2]])) ), 1)
    if (fc < 1) { fc <- paste("1/",fc,sep="") }
    test.stat <- wilcox.test(data.tmp[group==groupnames[1] ], data.tmp[group==groupnames[2] ], exact=F)$statistic
    c.index <- round(test.stat/(n.group[1]*n.group[2]),3)
    
    text.line <- NULL
    if(print.fc)      text.line <- paste(text.line, "  FC=", fc, sep="")
    if(print.average) text.line <- paste(text.line, "  mean=", round(average,2), sep="")
    if(print.sd)      text.line <- paste(text.line, "  SD=", round(sdev, 2), sep="")
    if(print.qvalue)  text.line <- paste(text.line, "  q-value=", round(qvalue, 4), sep="")
    if(print.cindex)  text.line <- paste(text.line, "  c-index=", c.index, sep="")
    title(sub=text.line, font=1, cex=1)
  }
}

plot.profile <- function(eset, group, geneid, legend.pos="topright")
{
  ####################################################################################
  ###                                                                              ###
  ### FUNCTION: Profile plot of a number of pre-specified genes                    ###
  ###                                                                              ###
  ### INPUT:                                                                       ###
  ###   eset          = matrix with gene expressions                               ###
  ###   group         = group membership                                           ###
  ###   genid         = which gene should be plotted                               ###
  ###   legend.pos    = position of legend in the plot                             ###        
  ###                                                                              ###
  ####################################################################################
  
  groupnames <- dimnames(table(group))
  groupnames <- c( groupnames[[1]][1], groupnames[[1]][2])
  n.group <- as.numeric(table(group))
  n.obs <- n.group[1]+n.group[2]
  data.tmp <- eset[geneid,]
  
  cols <- rainbow(length(geneid))
  
  plot(c(1, n.obs), c(floor(min(data.tmp)), ceiling(max(data.tmp))), col="white", bty="n", xlab="", ylab="", main="Profile plot", axes=FALSE, cex=2)
  axis(1, at=1:n.obs, labels=1:n.obs, cex=1)
  axis(1, at=c(n.group[1]/2,n.group[1]+(n.group[2]/2)), labels=groupnames, cex=1, line=2, tick=FALSE)
  axis(2, cex=1)
  box()
  abline(v=n.group[1]+0.5, col="gray", lwd=2)
  for (i in 1:length(geneid))
  {
    lines(1:n.obs, data.tmp[i,], col=cols[i])
    points(1:n.obs, data.tmp[i,], col=cols[i], pch=20)
  }
  legend(legend.pos, legend=geneid, lwd=1, col=cols, inset=0.02, lty=1, pch=20)
}

sensitivity.analyis <- function(input.folder, file.names, phenotype, n.g1, n.g2, seed1=712, seed2=524643181, IQR.cutoff=1)
{
  ####################################################################################  
  ###                                                                              ###
  ### FUNCTION: to compute a sensitivity analysis. One can change the arrays       ###
  ###           included in an analysis, or the IQR filter criterion               ###
  ###                                                                              ###
  ### At least one gene has to be selected with each method, otherwise you         ###
  ### will get some error.                                                         ###
  ###                                                                              ###
  ### INPUT:                                                                       ###
  ###   input.folder    = path where the arrays are stored                         ###
  ###   file.names      = file names of the arrays to be included                  ###
  ###   phenotype       = group membership for each array                          ###
  ###   n.g1            = number of samples in group 1                             ###
  ###   n.g2            = number of samples in group 2                             ###
  ###   seed1           = random number for SAM analysis                           ###
  ###   seed2           = random number for RankProduct-statistic                  ###
  ###   IQR.cutoff      = cut-off used for unspecific filtering with IQR           ###
  ###                                                                              ###
  ### OUTPUT: Matrix with an indicator for each gene stating, if the gene is       ###
  ###         selected with the filtering criterion ("filtered"), if the gene is   ###
  ###         selected as differentially expressed with SAM ("sam"), moderated     ###
  ###         t-statistic ("modt") and RankProduct-statistic ("RP").               ###
  ###         0 = not selected / 1 = selected                                      ###
  ###                                                                              ###
  ####################################################################################
  
  # Additionl R function: genelist.sam                                         
  genelist.sam2 <- function(x, y, fdr=0.01, nperms=100, t.sam=NA, delta.table=NA, refine.delta=TRUE,...)
  {
    ##############################################################################
    ###                                                                        ###
    ### FUNCTION: compute SAM analysis and save results                        ###
    ###                                                                        ###
    ### genelist.sam2 can be re-used with existing t.sam or delta.table        ###
    ### objects to save time                                                   ###
    ###                                                                        ###
    ### INPUT:                                                                 ###
    ###   x            = filtered gene expression data (exprs(filtered.data))  ###
    ###   y            = group membership (c(rep(1,nsBLC),rep(2,nsnBLC)))      ###
    ###   fdr          = favored FDR (default 0.01)                            ###
    ###   nperms       = number of permutations for SAM (default 100)          ###
    ###   t.sam        = object computed by samr function (optional)           ###
    ###   delta.table  = delta table (optional)                                ###
    ###   refine.delta = refine values for delta, TRUE or FALSE                ###
    ###                                                                        ###
    ### OUTPUT: a list with genelist, t.sam, delta.table and pi0.              ###
    ###                                                                        ###
    ##############################################################################
    
    data.sam <- list(x=x, geneid=as.character(1:nrow(x)),genenames=rownames(x), y=y,logged2=T)
    if(is.na(t.sam) | is.na(delta.table))
    {
      t.sam <- samr(data.sam, resp.type="Two class unpaired")
      if(!refine.delta) nvals=50
      else nvals=30
      delta.table <- samr.compute.delta.table(t.sam, nvals=nvals)      # use rough delta table to begin
    }
    delta.max <- delta.table[delta.table[,6] < fdr,1][1]               # select max delta for which FDR is controlled (90th percentile of FDR < FDRthreshold)
    delta.min <- delta.table[delta.table[,5] > 2*fdr, 1]               # select min delta for which FDR is controlled (median FDR > 2* FDRthreshold)
    delta.min <- as.vector(na.omit(delta.min))
    delta.min <- delta.min[length(delta.min)]
    # compute delta table once again, with more refined values for delta
    if(refine.delta == FALSE)
    {
      DEG.SAM <- samr.compute.siggenes.table(t.sam, delta.max, data.sam, delta.table, compute.localfdr=TRUE)
    }
    else
    {
      length.delta.table <- sum(-1*(is.na(delta.table[, 5]))+1)
      dels <- delta.min+(delta.max-delta.min)/length.delta.table*(0:length.delta.table)
      
      delta.table.fine <- samr.compute.delta.table(t.sam, dels=dels)    # use rough delta table to begin
      delta <- delta.table.fine[delta.table.fine[,5] < fdr,1][1]        # select delta for which FDR is controlled (median percentile of FDR < 0.01)
      DEG.SAM <- samr.compute.siggenes.table(t.sam, delta, data.sam, delta.table.fine, compute.localfdr=TRUE)
    }
    
    genelist <- rbind(DEG.SAM$genes.up, DEG.SAM$genes.lo)
    genelist[rank(-abs(as.numeric(genelist[,4]))),] <- genelist
    ngenes <- nrow(genelist)
    cat("Number of genes entered in analysis: ", nrow(x), "\n")
    cat("Estimated proportion of non-differentially expressed genes: ", t.sam$pi0, ".\n")
    cat("SAM found ", ngenes, " differentially expressed genes at an FDR of ", fdr, ".\n")
    list(genelist=genelist, t.sam=t.sam, delta.table=delta.table, pi0=t.sam$pi0)
  }
  
  #----------------------------------------------------------------------------------#
  library(affy)                              # read CEL files into R & preprocessing
  library(genefilter)                        # non-specific filtering
  library(limma)                             # gene expression analysis
  library(samr)                              # SAM
  library(RankProd)                          # RandProduct statistic
  
  memory.limit(2048)
  
  #----------------------------------------------------------------------------------#
  #    Data preprocessing                                                            #
  #----------------------------------------------------------------------------------#
  # Read data in from .cel files
  file.names2 <- paste(input.folder, file.names, sep="/")
  raw.data <- ReadAffy(filenames=file.names2)                             # this may take a few minutes
  
  tumor <- data.frame(SampleID=file.names, Tumor=phenotype)               # matrix indicating which kind of tumnor each file contains
  rownames(tumor) <- tumor$SampleID
  m <- match(rownames(tumor),sampleNames(raw.data))
  v <- data.frame(labelDescription=c("SampleID","Tumor"))
  phenoData(raw.data) <- new("AnnotatedDataFrame", data=data.frame(tumor[m,]), varMetadata=v)
  
  # Further standard settings
  raw.data2 <- exprs(raw.data)                                            # extracts expression data
  groupnames <- names(table(phenotype))
  
  #----------------------------------------------------------------------------------#
  #    Data normalization using RMA (background correction, normalization and        #
  #    calculation of summarized log-2 expression)                                   #
  #----------------------------------------------------------------------------------#
  rma.data <- rma(raw.data)
  rma.data2 <- exprs(rma.data)
  
  #----------------------------------------------------------------------------------#
  #    Unspecific filtering (based on IQR)                                           #
  #----------------------------------------------------------------------------------#
  filtered.data <- varFilter(rma.data, var.func=IQR, var.cutoff=IQR.cutoff, filterByQuantile=FALSE)
  eset <- exprs(filtered.data)
  
  #----------------------------------------------------------------------------------#
  #    Compute test statistics                                                       #
  #----------------------------------------------------------------------------------#
  # 1.) Moderated t-statistic (limma, another parametrization for the linear model)
  design <- model.matrix(~ 0+filtered.data$Tumor)
  colnames(design) <- groupnames
  fit <-lmFit(filtered.data, design)
  cont.matrix <- makeContrasts(contrasts=paste(groupnames[1], groupnames[2], sep="-"), levels=design)
  fit2 <- eBayes(contrasts.fit(fit, cont.matrix))
  x <- topTable(fit2, coef=1, adjust="BH", p.value=0.01, number=3000)
  
  # 2.) SAM
  set.seed(seed1)
  deg.sam <- genelist.sam2(x=exprs(filtered.data), y=c(rep(1,n.g1), rep(2,n.g2)), fdr=0.01)
  
  # 3.) RankProduct
  RP.out <- RP(exprs(filtered.data), cl=(phenotype==groupnames[2])+0, rand=seed2)
  RP.table <- topGene(RP.out, cutoff=0.005, method="pfp", logged=TRUE, logbase=2, gene.names=rownames(exprs(filtered.data)))
  RP.genelist <- rbind(RP.table$Table1, RP.table$Table2)
  
  #----------------------------------------------------------------------------------#
  #    Compare results of different test statistics                                  #
  #----------------------------------------------------------------------------------#
  # Find out which genes are contained in all lists
  sam.list <- deg.sam$genelist[,2]
  modt.list <- x[,1]
  RP.list <- rownames(RP.genelist)
  
  RPsam <- RP.list[!is.na(match(RP.list, sam.list))]
  
  #----------------------------------------------------------------------------------#
  #    Output                                                                        #
  #----------------------------------------------------------------------------------#
  results.sensi <- matrix(0, nrow=nrow(rma.data2), ncol=4,
                          dimnames=list(dimnames(rma.data2)[[1]], c("filtered", "sam", "modt", "RP")))
  results.sensi[match(dimnames(eset)[[1]], dimnames(results.sensi)[[1]]), "filtered"] <- 1
  results.sensi[match(sam.list, dimnames(results.sensi)[[1]]), "sam"] <- 1
  results.sensi[match(modt.list, dimnames(results.sensi)[[1]]), "modt"] <- 1
  results.sensi[match(RP.list, dimnames(results.sensi)[[1]]), "RP"] <- 1
  
  paste("From", nrow(rma.data2), "genes" , nrow(eset), "pass the IQR-filter criterion.", sep=" ")
  
  cat("\nNumber of genes selected by SAM and RankProduct-statistic:\n")
  print(sum(!is.na(match(RP.list, sam.list))))
  cat("\nNumber of genes selected by moderated t-statistic and RankProduct-statistic:\n")
  print(sum(!is.na(match(RP.list, modt.list))))
  cat("\nNumber of genes selected by SAM and moderated t-statistic\n")
  print(sum(!is.na(match(modt.list, sam.list))))
  cat("\nNumber of genes selected by SAM, moderated t-statistic and RankProduct-statistic:\n")
  print(sum(!is.na(match(RPsam, modt.list))))
  
  results.sensi
}

#-----------------------------------------------------------------------------------------------------------#
#    Data preprocessing                                                                                     #
#-----------------------------------------------------------------------------------------------------------#
# Read data from CEL files
file.names1 <- dir(input.folder)
file.names2 <- paste(input.folder,file.names1,sep="/")
raw.data <- ReadAffy(filenames=file.names2)                      # this may take a few minutes

nsBLC  <- as.numeric(table(phenotype)["BLC"])                    # number samples BLC group
nsnBLC <- as.numeric(table(phenotype)["nBLC"])                   # number samples nonBLC group

tumor <- data.frame(SampleID=file.names1, Tumor=phenotype)       # matrix indicating which kind of tumor each file contains
rownames(tumor) <- tumor$SampleID
m <- match(rownames(tumor),sampleNames(raw.data))
v <- data.frame(labelDescription=c("SampleID","Tumor"))
phenoData(raw.data) <- new("AnnotatedDataFrame", data=data.frame(tumor[m,]), varMetadata=v)
table(pData(raw.data)$Tumor)

# Further standard settings
memory.limit(2048)                                               # increases the memory limit

raw.data2 <- exprs(raw.data)                                     # extracts expression data
groupnames <- names(table(phenotype))
cols <- brewer.pal(4, "Dark2")                                   # nice colors for plots
rm(m, v, tumor)


#-----------------------------------------------------------------------------------------------------------#
#    Quality checks before normalization                                                                    #
#-----------------------------------------------------------------------------------------------------------#
# Plot of individual arrays on different scales (exemplified with the first array)
# The plot can be used to detect spatial artifacts. Of course, this should be done with all arrays.
palette.gray <- c(rep(gray(0:10/10), times=seq(1,41, by=4)))

array.no <- 1                                                           # request image plot of array number 1
jpeg(file=paste(output.folder, "1a - Quality before norm - Image of probe level data (first array).jpeg", sep="/"), quality=75)
image(raw.data[,array.no], transfo=function(x) x, col=palette.gray)   # natural scale
dev.off()
jpeg(file=paste(output.folder, "1b - Quality before norm - Image of probe level data (first array).jpeg", sep="/"), quality=75)
image(raw.data[,array.no], col=palette.gray)                          # log-2 transformed intensities
dev.off()
rm(palette.gray, array.no)

# Boxplot of probe level data for each array
jpeg(file=paste(output.folder, "2 - Quality before norm - Boxplots of probe level data.jpeg", sep="/"), quality=75)
boxplot(raw.data, main="Probe level data", col=c(rep(cols[1], 8), rep(cols[2], 8)), cex.axis=0.75)
axis(side=1, line=2, at=c(nsBLC/2, nsBLC+(nsnBLC/2)), labels=c(groupnames[1], groupnames[2]), tick=FALSE, font=2, cex=1.2)
dev.off()

# Plot of density estimates of probe level data for each array
jpeg(file=paste(output.folder, "3 - Quality before norm - Density estimates of probe level data.jpeg", sep="/"), quality=75)
hist(raw.data, col=cols, xlab="log2 intensities")
dev.off()

# MA-Plot to compare the probe level data of different arrays. We can either look at pairwise comparisons (pairs=TRUE)
# or comparisons to a reference array (default). We can either plot every point (default), or use a smoothing
# function (plot.method="smoothScatter").
jpeg(file=paste(output.folder, "4a - Quality before norm - MAplot of probe level data.jpeg", sep="/"), width=960, height=960, quality=75)
par(mfrow=c(4,4))                                                              
MAplot(raw.data, cex=0.75, cex.main=0.75)
mtext("M", 2, outer=TRUE)
mtext("A", 2, outer=TRUE)
dev.off()
jpeg(file=paste(output.folder, "4b - Quality before norm - MAplot of probe level data.jpeg", sep="/"), width=960, height=960, quality=75)
par(mfrow=c(4,4))
MAplot(raw.data, cex=0.75, plot.method="smoothScatter", cex.main=0.75)
mtext("M", 2, outer=TRUE)
mtext("A", 2, outer=TRUE)
par(mfrow=c(1,1))
dev.off()

jpeg(file=paste(output.folder, "4c - Quality before norm - MAplot of probe level data (array 1 and 2).jpeg", sep="/"), width=960, height=960, quality=75)
MAplot(raw.data[,1:2], pairs=TRUE, plot.method="smoothScatter")
dev.off()

# Assessment of specific probes
# probeNames(raw.data)                                                    # get probe names
probe <- "1007_s_at"
probeintensities <- pm(raw.data, probe)

jpeg(file=paste(output.folder, "5a - Quality before norm - PM probe intensities.jpeg", sep="/"), quality=75)
matplot(probeintensities, type="l", xlab="probe no.", ylab="pm probe intensity", main=probe)
dev.off()

jpeg(file=paste(output.folder, "5b - Quality before norm - PM probe intensities.jpeg", sep="/"), quality=75)
matplot(t(probeintensities), type="l", xlab="array no.", ylab="pm probe intensity", main=probe)
dev.off()
rm(probe, probeintensities)

# Tip: For more information of quality checks look at the affy vignettes.

# Affymetrix quality assessment (quality measures should be similar across all arrays)                                                    
quality <- qc(raw.data)
avbg(quality)                                                             # average background 
sfs(quality)                                                              # scale factor
percent.present(quality)                     
ratios(quality)[,1:2]             

# Compute and view concordance correlation and Spearman correlation
# (especially the computation of the concordance correlation takes time)
spearman.cor <- cor(raw.data2,  method=c("spearman"))
round(spearman.cor, 3)
mean(spearman.cor)                                                        # mean correlation over all arrays
round(apply(spearman.cor, 2, median), 3)                                  # median correlation per array

jpeg(file=paste(output.folder, "6 - Quality before norm - Spearman correlation of probe level data.jpeg", sep="/"), quality=75)
heatmap(x=-spearman.cor, main="Spearman correlation of probe level data", scale="none", col=heat.colors(256))
dev.off()

data.tmp1 <- spearman.cor
diag(data.tmp1) <- NA
jpeg(file=paste(output.folder, "7 - Quality before norm - Spearman correlation of probe level data.jpeg", sep="/"), quality=75)
plot.cor.boxplot(plot.data=data.tmp1, plot.title="Spearman correlation of probe level data", plot.cols=cols, 
                 n.g1=nsBLC, n.g2=nsnBLC, plot.groupnames=groupnames)
dev.off()
rm(data.tmp1)

concordance.cor <- matrix(NA, nrow=nsnBLC+nsBLC, ncol=nsnBLC+nsBLC, dimnames=list(dimnames(raw.data2)[[2]], dimnames(raw.data2)[[2]]))
for (i1 in 1:(nsnBLC+nsBLC-1)) {
  for (i2 in (i1+1):(nsnBLC+nsBLC)) {
    data.tmp2 <- na.omit(raw.data2[,c(i1,i2)])
    concordance.cor[i2, i1] <-
      concordance.cor[i1, i2] <- as.numeric(epi.ccc(x=data.tmp2[,1], y=data.tmp2[,2], ci="z-transform", conf.level=0.95)$rho.c[1])
  }
}
rm(data.tmp2)

diag(concordance.cor) <- 1
round(concordance.cor, 3)
mean(concordance.cor)
round(apply(concordance.cor, 2, median), 3)

jpeg(file=paste(output.folder, "8 - Quality before norm - Concordance correlation of probe level data.jpeg", sep="/"), quality=75)
heatmap(x=-concordance.cor, main="Concordance correlation of probe level data", scale="none", col=heat.colors(256))
dev.off()

data.tmp3 <- concordance.cor
diag(data.tmp3) <- NA
jpeg(file=paste(output.folder, "9 - Quality before norm - Concordance correlation of probe level data.jpeg", sep="/"), quality=75)
plot.cor.boxplot(plot.data=data.tmp3, plot.title="Concordance correlation of probe level data", plot.cols=cols,
                 n.g1=nsBLC, n.g2=nsnBLC, plot.groupnames=groupnames)
dev.off()
rm(data.tmp3)

# Fit robust linear model
dataPLM <- fitPLM(raw.data)
jpeg(file=paste(output.folder, "10 - Quality before norm - Robust linear model.jpeg", sep="/"), width=960, height=700, quality=75)
par(mfrow=c(2,1))
boxplot(dataPLM,main="NUSE")
Mbox(dataPLM,main="RLE")
dev.off()
par(mfrow=c(1,1))
rm(dataPLM)

# Compare two specific arrays: E. g. GSM85499.CEL and GSM85500.CEL
data.tmp4 <- na.omit(raw.data2[,c(13, 14)])
tmp.ccc <- epi.ccc(x=data.tmp4[,1], y=data.tmp4[,2], ci="z-transform", conf.level=0.95)
jpeg(file=paste(output.folder, "11 - Quality before norm - Concordance correlation of probe level data.jpeg", sep="/"), quality=75)
lab <- paste("CCC: ", round(tmp.ccc$rho.c[,1], digits = 2), " (95% CI ",
             round(tmp.ccc$rho.c[,2], digits = 2), " - ", round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = "")
par(pty = "s")
plot(data.tmp4[,1], data.tmp4[,2], xlab=dimnames(data.tmp4)[[2]][1], ylab=dimnames(data.tmp4)[[2]][2],
     pch=".", col="darkgray", cex=4)
abline(a=0, b=1, lty=1, col="black", lwd=2)
legend("topleft", inset=0.01, legend = c("Line of perfect concordance"),
       lty=1, lwd=2, bty="n", title=lab)
dev.off()
rm(data.tmp4, tmp.ccc, lab)


#-----------------------------------------------------------------------------------------------------------#
#    Data normalization using RMA (background correction, normalization and                                 #
#    calculation of summarized log-2 expression)                                                            #
#-----------------------------------------------------------------------------------------------------------#
rma.data  <- rma(raw.data)
rma.data2 <- exprs(rma.data)                                              # extracts gene expression data


#-----------------------------------------------------------------------------------------------------------#
#    Quality checks after normalization                                                                     #
#-----------------------------------------------------------------------------------------------------------#
# Boxplot
jpeg(file=paste(output.folder, "12 - Quality after norm - Boxplots of expression data.jpeg", sep="/"), quality=75)
boxplot(rma.data, main="Probe level data", col=c(rep(cols[1], 8), rep(cols[2], 8)), cex.axis=0.75)
axis(side=1, line=2, at=c(nsBLC/2, nsBLC+(nsnBLC/2)), labels=c(groupnames[1], groupnames[2]), tick=FALSE, font=2, cex=1.2)
dev.off()

# Plot of density estimates
jpeg(file=paste(output.folder, "13 - Quality after norm - Density estimates of expression data.jpeg", sep="/"), quality=75)
hist(rma.data, col=cols, xlab="log2 intensities")
dev.off()

# MA-Plot (examplified with 2/4 arrays)
jpeg(file=paste(output.folder, "14a - Quality after norm - MAplot of expression data.jpeg", sep="/"), width=700, height=700, quality=75)
par(mfrow=c(2,2))
MAplot(rma.data, cex=0.75, cex.main=0.75)
mtext("M", 2, outer=TRUE)
mtext("A", 2, outer=TRUE)
dev.off()
jpeg(file=paste(output.folder, "14b - Quality after norm - MAplot of expression data.jpeg", sep="/"), quality=75)
par(mfrow=c(2,2))
MAplot(rma.data, cex=0.75, plot.method="smoothScatter", cex.main=0.75)
mtext("M", 2, outer=TRUE)
mtext("A", 2, outer=TRUE)
dev.off()
par(mfrow=c(1,1))

jpeg(file=paste(output.folder, "14c - Quality after norm - MAplot of expression data.jpeg", sep="/"), quality=75)
MAplot(rma.data[,1:2], pairs=TRUE, plot.method="smoothScatter")         # MAplot of the first two arrays
dev.off()
#jpeg(file=paste(output.folder, "14d - Quality after norm - MAplot of expression data.jpeg", sep="/"), quality=75)
#  MAplot(rma.data, cex=0.75, plot.method="smoothScatter", cex.main=0.75, pairs=TRUE)
#dev.off()

# Compute and view concordance correlation and Spearman correlation
# (especially the computation of the concordance correlation takes time)
spearman.cor <- cor(rma.data2,  method=c("spearman"))
round(spearman.cor, 3)
mean(spearman.cor)                                                        # mean correlation over all arrays
round(apply(spearman.cor, 2, median), 3)                                  # median correlation per array

jpeg(file=paste(output.folder, "15 - Quality after norm - Spearman correlation of expression data.jpeg", sep="/"), quality=75)
heatmap(x=-spearman.cor, main="Spearman correlation of expression data", scale="none", col=heat.colors(256))
dev.off()

data.tmp1 <- spearman.cor
diag(data.tmp1) <- NA
jpeg(file=paste(output.folder, "16 - Quality after norm - Spearman correlation of expression data.jpeg", sep="/"), quality=75)
plot.cor.boxplot(plot.data=data.tmp1, plot.title="Spearman correlation of expression data", plot.cols=cols,
                 n.g1=nsBLC, n.g2=nsnBLC, plot.groupnames=groupnames)
dev.off()
rm(data.tmp1)

concordance.cor <- matrix(NA, nrow=nsnBLC+nsBLC, ncol=nsnBLC+nsBLC, dimnames=list(dimnames(rma.data2)[[2]], dimnames(rma.data2)[[2]]))
for (i1 in 1:(nsnBLC+nsBLC-1)) {
  for (i2 in (i1+1):(nsnBLC+nsBLC)) {
    data.tmp2 <- na.omit(rma.data2[,c(i1,i2)])
    concordance.cor[i2, i1] <-
      concordance.cor[i1, i2] <- as.numeric(epi.ccc(x=data.tmp2[,1], y=data.tmp2[,2], ci="z-transform", conf.level=0.95)$rho.c[1])
  }
}
rm(data.tmp2, i1, i2)

diag(concordance.cor) <- 1
round(concordance.cor, 3)
mean(concordance.cor)
round(apply(concordance.cor, 2, median), 3)

jpeg(file=paste(output.folder, "17 - Quality after norm - Concordance correlation of expression data.jpeg", sep="/"), quality=75)
heatmap(x=-concordance.cor, main="Concordance correlation of expression data", scale="none", col=heat.colors(256))
dev.off()

data.tmp3 <- concordance.cor
diag(data.tmp3) <- NA
jpeg(file=paste(output.folder, "18 - Quality after norm - Concordance correlation of expression data.jpeg", sep="/"), quality=75)
plot.cor.boxplot(plot.data=data.tmp3, plot.title="Concordance correlation of expression data", plot.cols=cols,
                 n.g1=nsBLC, n.g2=nsnBLC, plot.groupnames=groupnames)
dev.off()
rm(data.tmp3)

# Spectral map of first two principal components
sma <- mpm(data=data.frame(rownames(rma.data2), rma.data2), logtrans=FALSE, center="row",
           normal="none")

jpeg(file=paste(output.folder, "19 - Quality after norm - Spectral map.jpeg", sep="/"), quality=75)
r <- plot(sma, label.tol=0, scale="uvc", dim=c(1,2), col.size=2, col.areas=FALSE,
          zoom = c(1,1.2), col.group=(phenotype=="BLC")+1)
dev.off()
rm(r, sma)

# Compare two specific arrays: E. g. GSM85499.CEL and GSM85500.CEL
data.tmp4 <- na.omit(rma.data2[,c(13, 14)])
tmp.ccc <- epi.ccc(x=data.tmp4[,1], y=data.tmp4[,2], ci="z-transform", conf.level=0.95)
jpeg(file=paste(output.folder, "20 - Quality after norm - Concordance correlation of expression data.jpeg", sep="/"), quality=75)
lab <- paste("CCC: ", round(tmp.ccc$rho.c[,1], digits = 2), " (95% CI ",
             round(tmp.ccc$rho.c[,2], digits = 2), " - ", round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = "")
par(pty = "s")
plot(data.tmp4[,1], data.tmp4[,2], xlab=dimnames(data.tmp4)[[2]][1], ylab=dimnames(data.tmp4)[[2]][2],
     pch=".", col="darkgray", cex=4)
abline(a=0, b=1, lty=1, col="black", lwd=2)
legend("topleft", inset=0.01, legend = c("Line of perfect concordance"),
       lty=1, lwd=2, bty="n", title=lab)
dev.off()
rm(data.tmp4, tmp.ccc, lab)

#-----------------------------------------------------------------------------------------------------------#
#    Discard arrays (based on the results of the quality checks)                                            #
#-----------------------------------------------------------------------------------------------------------#
# Discard arrays with low correlation to other arrays
round(sort(apply(spearman.cor, 1, min)), digits=5)                       # gives the minimum spearman correlation per array
round(sort(apply(spearman.cor, 1, median)), digits=5)                    # gives the median spearman correlation per array   

round(sort(apply(concordance.cor,1, min)), digits=5)
round(sort(apply(concordance.cor,1, median)), digits=5)


#-----------------------------------------------------------------------------------------------------------#
#    Unspecific filtering (based on IQR)                                                                    #
#-----------------------------------------------------------------------------------------------------------#
variation <- apply(rma.data2, 1, IQR)
jpeg(file=paste(output.folder, "21 - Variation of expression data.jpeg", sep="/"), quality=75)
tmp <- summary(variation)
boxplot(variation, main="IQR of expression data", sub=paste(names(tmp)[1],
                                                            as.numeric(round(tmp[1], 3)), names(tmp)[2], as.numeric(round(tmp[2], 3)),
                                                            names(tmp)[3], as.numeric(round(tmp[3], 3)), names(tmp)[4], as.numeric(round(tmp[4], 3)),
                                                            names(tmp)[5], as.numeric(round(tmp[5], 3)), names(tmp)[6], as.numeric(round(tmp[6], 3)), sep=" "))
dev.off()
rm(tmp)

# Filtering: var.cutoff is the filtering cut-off
filtered.data <- varFilter(rma.data, var.func=IQR, var.cutoff=1, filterByQuantile=FALSE)
# reduces number of genes to 5874 - only these have an IQR>1!
eset <- exprs(filtered.data)
paste("From", nrow(rma.data2), "genes" , nrow(eset), "pass the IQR-filter criterion.", sep=" ")


#-----------------------------------------------------------------------------------------------------------#
#    Visualizating differential expression                                                                  #
#-----------------------------------------------------------------------------------------------------------#
# Principal components and Spectral map
sma <- mpm(data=data.frame(rownames(eset), eset), logtrans=FALSE, center="row", normal="none")
sma2 <- summary(sma)

# Figure 1 from the book
jpeg(file=paste(output.folder, "22 - Plot of first 2 principal components.jpeg", sep="/"), quality=75)
plot(sma2$Columns[,3], sma2$Columns[,4], pch=(phenotype=="BLC")+1, xlab=paste("PC1 (", round(sma2$VPF[1,1]*100,0),"%)"),
     ylab=paste("PC2 (",round(sma2$VPF[1,2]*100,0),"%)"))
dev.off()

jpeg(file=paste(output.folder, "23a - Spectral map.jpeg", sep="/"), quality=75)
r <- plot(sma, label.tol=0, scale="uvc", dim=c(1,2), col.size=2, col.areas=FALSE,
          zoom = c(1,1.2), col.group=(phenotype=="BLC")+1)
dev.off()
jpeg(file=paste(output.folder, "23b - Spectral map.jpeg", sep="/"), quality=75)
r <- plot(sma, label.tol=0, scale="uvc", dim=c(1,2), col.size=2, col.areas=FALSE,
          col.group=(phenotype=="BLC")+1, zoom = c(1,1.2), do.smoothScatter=TRUE)
dev.off()                                                                 # Spectral map of first two principal components

jpeg(file=paste(output.folder, "23c - Spectral map.jpeg", sep="/"), quality=75)
r <- plot(sma, label.tol=25, scale = "uvc", dim=c(1,2), col.size=2, col.areas=FALSE,
          lab.size=0.7, col.group=(phenotype=="BLC")+1)
dev.off()

# Give the 25 most extreme genes, i. e. most distant from center
subset(as.data.frame(r$Rows), Select==1)
rm(r, sma, sma2)

#-----------------------------------------------------------------------------------------------------------#
#    Statistics comparison: compute test statistics, correct p-values and                                   #
#    store lists of differentially expressed genes                                                          #
#-----------------------------------------------------------------------------------------------------------#
# 1.) moderated t-statistic (limma, another parametrization for the linear model)
design <- model.matrix(~ 0+filtered.data$Tumor)
colnames(design) <- groupnames
fit <-lmFit(filtered.data, design)
cont.matrix <- makeContrasts(BLC-nBLC, levels=design)
fit2 <- eBayes(contrasts.fit(fit, cont.matrix))
x <- topTable(fit2, coef=1, adjust="BH", p.value=0.01, number=3000)
nrow(x)
write.table(data.frame(x), paste(output.folder, "DEG_limma.txt", sep="/"), quote=F,sep="\t", row.names=F, dec=",")
xfull <- topTable(fit2, coef=1, adjust="BH", p.value=1, number=6000)

# 2.) SAM
set.seed(7123981)
deg.sam <- genelist.sam(x=exprs(filtered.data), y=c(rep(1,nsBLC), rep(2,nsnBLC)), fdr=0.01)

cat("\nSAM: Delta Table:\n")
deg.sam$delta.table
cat("\nSAM: Delta for which FDR of 0.01 is controlled (90th percentile of FDR < 0.01):\n")
deg.sam$delta

jpeg(file=paste(output.folder, "24 - SAM plot.jpeg", sep="/"))
samr.plot(deg.sam$t.sam, deg.sam$delta)
dev.off()

write.table(data.frame(deg.sam$genelist), paste(output.folder, "DEG_SAM.txt", sep="/"),
            quote=F, sep="\t", row.names=F, dec=",")

siggenes.table <- samr.compute.siggenes.table(samr.obj=deg.sam$t.sam, del=deg.sam$delta, data=deg.sam$data.sam, 
                                              delta.table=deg.sam$delta.table, min.foldchange=3, all.genes=FALSE, compute.localfdr=TRUE)
# Here we use a large min.foldchange in order to get only a short list of genes.
cat("\nSAM: Significant genes\n")
siggenes.table

cat("\nSAM: Miss rate\n")                               # shows the local false negative rate
samr.missrate(samr.obj=deg.sam$t.sam, del=deg.sam$delta, delta.table=deg.sam$delta.table)

# 3.) RankProduct
RP.out <- RP(exprs(filtered.data), cl=(phenotype=="nBLC")+0, rand=7123981)

# Identify differentially expressed genes
RP.table <- topGene(RP.out, cutoff=0.005, method="pfp", logged=TRUE, logbase=2, gene.names=rownames(exprs(filtered.data)))
RP.genelist <- rbind(RP.table$Table1, RP.table$Table2)
cat("\nRankProduct: number of differentially expressed genes\n")
nrow(RP.genelist)
write.table(data.frame(RP.genelist), paste(output.folder, "DEG_RP.txt", sep="/"), quote=F, sep="\t", row.names=T, dec=",")

#-----------------------------------------------------------------------------------------------------------#
#    Visualize results                                                                                      #
#-----------------------------------------------------------------------------------------------------------#
# Cluster analysis: Heat map of the first 100 genes selected as differentially expressed by limma.
DEG.genes <- match(x$ID[1:100], dimnames(eset)[[1]])
data.tmp6 <- eset[DEG.genes, ]
jpeg(file=paste(output.folder, "25a - Heatmap (limma results).jpeg", sep="/"), quality=75)
heatmap(x=data.tmp6, xlab=" ", ylab="gene expressions", main="Differentially expressed genes", scale="none", labRow=NA,
        cexCol=0.75, margins=c(5, 2))
dev.off()

jpeg(file=paste(output.folder, "25b - Heatmap (limma results).jpeg", sep="/"), quality=75)
heatmap(x=data.tmp6, Rowv=NA, scale="none", main="Differentially expressed genes",
        labRow=NA, cexCol=0.75, margins=c(5, 2), xlab=" ", ylab="gene expressions")
dev.off()

# Look at specific genes: Dotplot
jpeg(file=paste(output.folder, "26a - Geneplot specific genes.jpeg", sep="/"), quality=75)
plot.genes(eset=eset, group=phenotype, geneid="205044_at", qtable=xfull, print.fc=TRUE, print.average=TRUE, print.sd=TRUE, print.cindex=TRUE, print.qvalue=TRUE)
dev.off()

jpeg(file=paste(output.folder, "26b - Geneplot specific genes.jpeg", sep="/"), quality=75)
par(mfrow=c(2,2), oma=c(0,0,0,0))
plot.genes(eset=eset, group=phenotype, geneid="205044_at",   qtable=xfull, print.cindex=TRUE)
plot.genes(eset=eset, group=phenotype, geneid="226684_at",   qtable=xfull, print.cindex=TRUE)
plot.genes(eset=eset, group=phenotype, geneid="220892_s_at", qtable=xfull, print.cindex=TRUE)
plot.genes(eset=eset, group=phenotype, geneid="219654_at",   qtable=xfull, print.cindex=TRUE)
dev.off()

jpeg(file=paste(output.folder, "26c - Geneplot specific genes.jpeg", sep="/"), quality=75)
plot.genes(eset=eset, group=phenotype, geneid="1553613_s_at")
plot.genes(eset=eset, group=phenotype, geneid="203963_at")
plot.genes(eset=eset, group=phenotype, geneid="215867_x_at")
plot.genes(eset=eset, group=phenotype, geneid="220425_x_at")
dev.off()

jpeg(file=paste(output.folder, "26d - Geneplot specific genes.jpeg", sep="/"), quality=75)
par(mfrow=c(1,2))
plot.genes(eset=eset, group=phenotype, geneid="224590_at")
plot.genes(eset=eset, group=phenotype, geneid="229927_at")
dev.off()
par(mfrow=c(1,1))

# Look at specific genes: Profile plot
jpeg(file=paste(output.folder, "27 – Profileplot.jpeg", sep="/"))
plot.profile(eset=eset, group=phenotype, geneid=c("205044_at", "1553613_s_at", "204667_at", "203963_at"))
dev.off()


#-----------------------------------------------------------------------------------------------------------#
#    Compare results of different test statistics                                                           #
#-----------------------------------------------------------------------------------------------------------#
# Find out which genes are contained in all lists
sam.list <- deg.sam$genelist[,2]
modt.list <- x[,1]
RP.list <- rownames(RP.genelist)
cat("\nNumber of genes selected by SAM and RankProduct-statistic:\n")
sum(!is.na(match(RP.list, sam.list)))
cat("\nNumber of genes selected by moderated t-statistic and RankProduct-statistic:\n")
sum(!is.na(match(RP.list, modt.list)))
cat("\nNumber of genes selected by SAM and moderated t-statistic\n")
sum(!is.na(match(modt.list, sam.list)))

# All three
RPsam <- RP.list[!is.na(match(RP.list, sam.list))]
cat("\nNumber of genes selected by SAM, moderated t-statistic and RankProduct-statistic:\n")
sum(!is.na(match(RPsam, modt.list)))

# We summarize some results: For each gene we store the following information in a matrix:
# "filtered"=1 == gene is filtered, 
# "tt"=1       == gene is selected by t-statistic,
# "modt"=1     == gene is selected by moderated t-statistic, 
# "RP"=1       == gene is selected by RankProduct-statistic
results.original <- matrix(0, nrow=nrow(rma.data2), ncol=4, dimnames=list(dimnames(rma.data2)[[1]], c("filtered", "sam", "modt", "RP")))
results.original[match(dimnames(eset)[[1]], dimnames(results.original)[[1]]), "filtered"] <- 1
results.original[match(sam.list,            dimnames(results.original)[[1]]), "sam"]  <- 1
results.original[match(modt.list,           dimnames(results.original)[[1]]), "modt"] <- 1
results.original[match(RP.list,             dimnames(results.original)[[1]]), "RP"]   <- 1

head(results.original)
write.table(data.frame(results.original), paste(output.folder, "results.original.txt", sep="/"), quote=F, sep="\t", row.names=T, dec=",")

#-----------------------------------------------------------------------------------------------------------#
#    Sensitivity analysis                                                                                   #
#-----------------------------------------------------------------------------------------------------------#
save.image(file=paste(output.folder, "example.RData", sep="/"))
keep(results.original, file.names1, input.folder, output.folder, sensitivity.analyis, sure=TRUE)

# No filter prior to selection
results.nofilter <- sensitivity.analyis(input.folder=input.folder, file.names=file.names1[], phenotype=c(rep("BLC", 8), rep("nBLC", 8)),
                                        n.g1=8, n.g2=8, seed1=6123981, seed2=5123981, IQR.cutoff=0)

cat("\nSAM: Number of DEG of original analysis vs. sensitivity analysis (no filter) (0=not selected 1=selected)\n")
table(results.original[,"sam"], results.nofilter[,"sam"], dnn=c("Original", "Sensitivity"))

cat("\nmodt: Number of DEG of original analysis vs. sensitivity analysis (no filter) (0=not selected 1=selected)\n")
table(results.original[,"modt"], results.nofilter[,"modt"], dnn=c("Original", "Sensitivity"))

cat("\nRP: Number of DEG of original analysis vs. sensitivity analysis (no filter) (0=not selected 1=selected)\n")
table(results.original[,"RP"], results.nofilter[,"RP"], dnn=c("Original", "Sensitivity"))

# Another cut-off for the IQR-filter-criterion: Here analysis with IQR cut-off of 1.2
results.IQRcut <- sensitivity.analyis(input.folder=input.folder, file.names=file.names1, phenotype=c(rep("BLC", 8), rep("nBLC", 8)),
                                      n.g1=8, n.g2=8, seed1=4123982, seed2=6523845, IQR.cutoff=1.2)

cat("\nSAM: Number of DEG of original analysis vs. sensitivity analysis (IQR cut-off 1.2) (0=not selected 1=selected)\n")
table(results.original[,"sam"], results.IQRcut[,"sam"], dnn=c("Original", "Sensitivity"))

cat("\nmodt: Number of DEG of original analysis vs. sensitivity analysis (IQR cut-off 1.2) (0=not selected 1=selected)\n")
table(results.original[,"modt"], results.IQRcut[,"modt"], dnn=c("Original", "Sensitivity"))

cat("\nRP: Number of DEG of original analysis vs. sensitivity analysis (IQR cut-off 1.2) (0=not selected 1=selected)\n")
table(results.original[,"RP"], results.IQRcut[,"RP"], dnn=c("Original", "Sensitivity"))

# Exclude specific arrays from the analysis: Here analysis without array ...
results.excludearray <- sensitivity.analyis(input.folder=input.folder, file.names=file.names1[-16], phenotype=c(rep("BLC", 8), rep("nBLC", 7)),
                                            n.g1=8, n.g2=7, seed1=3123983, seed2=58964752, IQR.cutoff=1)

cat("\nSAM: Number of DEG of original analysis vs. sensitivity analysis (array GSM85503.CEL is excluded) (0=not selected 1=selected)\n")
table(results.original[,"sam"], results.excludearray[,"sam"], dnn=c("Original", "Sensitivity"))

cat("\nmodt: Number of DEG of original analysis vs. sensitivity analysis (array GSM85503.CEL is excluded) (0=not selected 1=selected)\n")
table(results.original[,"modt"], results.excludearray[,"modt"], dnn=c("Original", "Sensitivity"))

cat("\nRP: Number of DEG of original analysis vs. sensitivity analysis (array GSM85503.CEL is excluded) (0=not selected 1=selected)\n")
table(results.original[,"RP"], results.excludearray[,"RP"], dnn=c("Original", "Sensitivity"))

#-----------------------------------------------------------------------------------------------------------#
#    Save workspace                                                                                         #
#-----------------------------------------------------------------------------------------------------------#
load(file=paste(output.folder, "example.RData", sep="/"))
save.image(file=paste(output.folder, "example.RData", sep="/"))

