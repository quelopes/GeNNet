enrichment_GOstats2 = function(listData,nameExp,annot.name,wayOut){
  
  print("*** ACTIVITY BIOLOGICAL PROCESS OVER REPRESENTED  ***")
  
  # wayOut = dirW
  # annot.name = orgAnnot
  # nameExp = nameExp

  mat.exp = listData[[6]]
  entrez =  listData[[5]]$entrez
  universe = names(listData[[2]][,1])
  
  
  # --- parameters ---
  hgCutoff = 0.001
  
  # ==============  
  # --- GoStat ---
  # ==============
  # cluster = 6 # = indice number of list. list[[6]]
  ncluster = max(mat.exp) # name of data in load 
  
  # check what col use (model or orthologo)
  check = universe[match(entrez,universe)]

  
  names(mat.exp) = entrez
  
  mat.GO = NULL  
  l.value = NULL
  
  for(k in 1:ncluster){   
    dados.sel = which(mat.exp == k)
    names.dados = rownames(mat.exp)
    names.dados.sel = names.dados[dados.sel]
    
    params = new("GOHyperGParams",
                 geneIds = names.dados.sel,
                 universeGeneIds = universe,
                 annotation = annot.name,
                 ontology="BP",
                 pvalueCutoff = hgCutoff,
                 conditional = FALSE, #TRUE
                 testDirection = "over")
    paramsCond = params
    conditional(paramsCond) = TRUE
    hgOver.BP = hyperGTest(params)
    
    if(nrow(summary(hgOver.BP)) >= 1){

      # Get the entrez id for each GO term:
      hgOver.BP.probes = geneIds(hgOver.BP)
      hgOver.BP.probes = geneIdsByCategory(hgOver.BP, catids = summary(hgOver.BP)$GOBPID)
      GOtable = data.frame(GOterm=summary(hgOver.BP)$GOBPID, term = summary(hgOver.BP)$Term,Pvalue=format(summary(hgOver.BP)$Pvalue, scientific=T), 
                           NoGOsize=paste(summary(hgOver.BP)$Count,summary(hgOver.BP)$Size, sep=";"))
      GOtable = cbind(GOtable,cluster = paste("C",k,nameExp,sep="_"))
      
      mat.agreg = NULL 
      
      for(j in 1:length(hgOver.BP.probes)){
        tmp = hgOver.BP.probes[[j]]
        e.go = entrez[which(is.na(match(entrez,tmp))==FALSE)]
        l.symbol = paste(e.go,";",sep="",collapse='')
        l.symbol = substr(l.symbol, 1, nchar(l.symbol)-1)
        
        mat.agreg = rbind(mat.agreg,l.symbol)
      }
      colnames(mat.agreg) = "entrez"
      GOtable = cbind(GOtable,mat.agreg)
      mat.GO = rbind(mat.GO,GOtable)
    }
  }
  
  mat.GO$term = gsub(",",";",mat.GO$term)
  write.table(mat.GO,file=paste(wayOut,nameExp,"-MatrixGO.csv",sep=""),sep=",", row.names=T, col.names=T,quote=F)
  
  names.list = c(names(listData),"functionalAnalises")
  listData[[7]] = mat.GO
  names(listData) = names.list
  
  listData
  
}
