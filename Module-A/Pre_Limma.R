pre_Limma = function(eSet,pValue,foldChange,wayOut,name.exp){
  
  print("*** ACTIVITY DIFFERENTIAL EXPRESSION ***")	
  
  # pValue = pValue
  # foldChange = FC
  # wayOut = dirW
  # name.exp = names(GeNNet[1])
  # eSet = GeNNet[[1]]
  
  # ---------------------	
  # *** MATRIX DESIGN *** 
  # ---------------------	
  targets = pData(eSet)
  
  f = targets$SETS
  f = factor(f)
  design = model.matrix(~-1+f)
  
  # f <- factor(targets$SETS)
  # design <- model.matrix(~f)
  # fit <- eBayes(lmFit(eSet,design))
  
  # *** NOMES UTILIZAVEIS COLUNAS ***
  namesTime = make.names(unique(f))
  
  colnames(design) = namesTime
  fit = lmFit(eSet, design)
  
  matExpProbe = fit$coefficients # valor ajustado, melhor que a media
  write.table(matExpProbe, paste(wayOut,name.exp,"-MatrixProbe-lmFit.csv",sep=""), sep=",", row.names=T, col.names=T, dec=".", quote=F)

  # *** MATRIX CONTRAST ***
  # diff = NULL 	# diff ser? o nome dado a cada contraste da matriz	
  # for(i in 1:(length(namesTime)-1)){
  #   diff = c(diff,paste(namesTime[i+1],namesTime[i],sep="-"))
  # }
  # cont.matrix = makeContrasts(contrasts=diff,levels=design)
  
  diff = NULL 	# diff ser? o nome dado a cada contraste da matriz	
  for(i in 1:(length(namesTime))){
    for(j in 2:(length(namesTime))){
      if(j > i){
        diff = c(diff,paste(namesTime[j],namesTime[i],sep="-"))
      }
    }
    rm(j)
  }
  rm(i)
  cont.matrix = makeContrasts(contrasts=diff,levels=design)
  
  
  # -------------------------------	
  # *** FITTING LINEAR FUNCTION *** 
  # -------------------------------	
  fit2  = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)
  
  results = decideTests(fit2)
  
  summary(results)
  result.fit = summary(results)
  #vennDiagram(results)
  
  lines.expr.set = nrow(exprs(eSet))
  
  # -------------------	
  # *** DIFFERENCES *** 
  # -------------------	
  
  for(k in 1:ncol(cont.matrix)){
    pdf(paste(wayOut,name.exp,"-VolcanoPlot-",diff[k],".pdf",sep="")) 
    
    diff.val = topTable(fit2, coef=diff[k],sort.by="none", number=lines.expr.set,adjust.method="BH")#, sort.by="none", adjust.method="BH", p.value=1, lfc=0)
    
    diff.val$odds = exp(diff.val[,"B"])
    diff.val$odds.Prob = (diff.val$odds)/(1+diff.val$odds)
    
    # Selecionando valores de acordo com o threshold
    ind = with(diff.val, which(adj.P.Val  <= pValue & abs(logFC) >= foldChange))
    diff.val.DE = diff.val[ind,]
    
    if(is.character(diff.val[1,1]) == TRUE){
      probes = as.character(diff.val.DE$ID)
      probes2 = as.character(diff.val$ID)
    }
    if(is.character(diff.val[1,1]) == FALSE){
      probes = rownames(diff.val.DE)
      probes2 = rownames(diff.val)
    }
    
    # Volcano
    volcanoplot(fit2, coef=diff[k], highlight=50, 
                # names=fit2$genes$ID, 
                main=paste("Volcano",diff[k],sep=": "))
    abline(h=-log10(0.05), lty=2)
    abline(v=-1.5, lty=2, col="green")
    abline(v=1.5, lty=2, col="red")
    
    if(k == 1){
      probes.unique = probes
      probes.num = length(probes)
      # probes.unique.all = all pvalue and FC 
      probes.unique.all = probes2
    }
    if(k > 1){
      probes.unique = unique(c(probes.unique,probes)) 
      probes.num = c(probes.num,length(probes))
      
      probes.unique.all = unique(c(probes.unique.all,probes2)) 
    } 
    dev.off()
  }
  
 # ---
  p.num = as.data.frame(probes.num)
  p.numNames = diff
  
  # rownames(p.num) = m2
  p.num = t(p.num)
  
  # *** MERGE PROBES COM FC AND PVALUE SELECTED ***
  vec = c(1:length(probes.unique))
  vec.all = c(1:length(probes.unique.all))
  
  
  if(length(probes.unique)!=0){
    # criando um data frame para armazenar valores de fold-change e de pvalue
    list.DE.all = data.frame(ID=probes.unique,vec)
    list.all = data.frame(ID=probes.unique.all,vec.all)
    
    for(m in 1:ncol(cont.matrix)){
      mat = topTable(fit2, coef=diff[m], number=lines.expr.set)
      if(is.character(diff.val[1,1]) == TRUE){
        list.DE.all = merge(list.DE.all, mat[,c(1,2,6)],by="ID",all.x=TRUE)
        mat = topTable(fit2, coef=diff[m], number=lines.expr.set)
        list.all = merge(list.all, mat[,c(1,2,6)],by="ID",all.x=TRUE)
      }
      if(is.character(diff.val[1,1]) == FALSE){
        list.DE.all = merge(list.DE.all, mat[,c(1,5)],by.x="ID",by.y='row.names',all.x=TRUE)
        mat = topTable(fit2, coef=diff[m], number=lines.expr.set)
        list.all = merge(list.all, mat[,c(1,5)],by.x="ID",by.y='row.names',all.x=TRUE)
      }
    }
    
    names = as.character(list.DE.all[,1])
    list.DE.all = list.DE.all[,-c(1,2)]
    rownames(list.DE.all) = names
    selecionados = rownames(list.DE.all)
  }
  
  
  if(length(probes.unique)==0){
    print("No probes/genes Differentially Expressed!")
    list.all = data.frame(ID=probes.unique.all,vec.all)
    
    for(m in 1:ncol(cont.matrix)){
      mat = topTable(fit2, coef=diff[m], number=lines.expr.set)
      if(is.character(diff.val[1,1]) == TRUE){
        mat = topTable(fit2, coef=diff[m], number=lines.expr.set)
        list.all = merge(list.all, mat[,c(1,2,6)],by="ID",all.x=TRUE)
      }
      if(is.character(diff.val[1,1]) == FALSE){
        mat = topTable(fit2, coef=diff[m], number=lines.expr.set)
        list.all = merge(list.all, mat[,c(1,5)],by.x="ID",by.y='row.names',all.x=TRUE)
      }
    }
    selecionados=NULL
  }
  
  names2 = as.character(list.all[,1])
  list.all = list.all[,-c(1,2)]
  rownames(list.all) = names2

  # ---------------------------------	
  # *** CHANGE IN EXPRS AND WRITE ***	
  # ---------------------------------
  full.list = list(eSet,matExpProbe,list.all,selecionados)
  names(full.list) = c(name.exp,"matExpProbe","mat.FC.pVal","selecionados")
  full.list
}
