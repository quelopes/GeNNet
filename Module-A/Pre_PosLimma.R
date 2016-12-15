pre_PosLimma = function(annotation,list.Exp,nameExp,dirOut){
  
  # annotation = mat_annot
  # list.Exp = GeNNet
  # names(list.Exp)
  # dim(mat_annot)
  # dirOut = dirW
  
  print("*** ACTIVITY ANNOTATION PROBES TO GENES ***")  
  
  expMat = list.Exp[[2]]
  nam.exp = names(list.Exp[1])
  
  # Order: Probe_Id [1], Entrez_Id [2], Gene_Symbol [3]
  probe = 1
  entrez = 2
  geneSymbol = 3
  
  # ==========
  # --- IF ---
  # ==========
  if("Probe.Set.ID" %in% colnames(annotation)){
    annotation = annotation[,c("Probe.Set.ID","Entrez.Gene","Gene.Symbol")]
    
    # REMOVE LINES WITH NULL
    
    sel = annotation[annotation$Entrez.Gene == "null" | annotation$Entrez.Gene == "---",] # add other characters <<<<<<<<<<<<<<<<<<<
    if(length(sel[,2]) > 0){
      dados.anot = annotation[-as.numeric(rownames(sel)),] 
    }
    
    # DOUBLE ENTREZ IN SOME LINE
    dados.anot[,entrez] = as.character(dados.anot[,entrez])
    dados.anot[,geneSymbol] = as.character(dados.anot[,geneSymbol])
    
    caracter1 = "-"
    caracter2 = "///"
    
    for(i in 1:nrow(dados.anot)){
      tm1 = length(grep(caracter1, dados.anot[i,entrez]))
      tm2 = length(grep(caracter2, dados.anot[i,entrez]))
      
      if(tm1 != 0 || tm2 != 0){
        if(tm1 != 0){
          get1 = strsplit(dados.anot[i,entrez], split=caracter1)[[1]][1]	
          get2 = strsplit(dados.anot[i,geneSymbol], split=caracter1)[[1]][1]	
          dados.anot[i,entrez] = get1
          dados.anot[i,geneSymbol] = get2
        }
        if(tm2 != 0){
          get1 = strsplit(dados.anot[i,entrez], split=caracter2)[[1]][1]
          get2 = strsplit(dados.anot[i,geneSymbol], split=caracter2)[[1]][1]	
          
          dados.anot[i,entrez] = gsub(" *$", "",get1)
          dados.anot[i,geneSymbol] = get2
        }
      }
    }
    expP = expMat[!is.na(match(row.names(expMat),dados.anot$Probe.Set.ID)),  ]
    rownames(expP)
    names = dados.anot[match(row.names(expP),dados.anot[,1]),"Entrez.Gene"]
    row.names(expP) = names
    # head(expP)
    # OBS: podemos observar aqui a duplicidade de entrezId nas linhas. Selecionei apenas o primeiro entrezId unico, mas futuramente considerar outros metodos.
    expP = expP[which(!duplicated(rownames(expP))),]
    write.table(expP, paste(dirOut,nam.exp,"-MatrixExpressionEntrez.csv",sep=""),sep=",", row.names=T, col.names=T, dec=".", quote=F)
    
    # Matriz diferencialmente expressa
    # exp.sel = matExpProbe[match(rownames(list.DE.all),row.names(matExpProbe)),  ]
    
    exp.sel = expMat[match(list.Exp[[4]],row.names(expMat)),  ]
    # head(exp.sel)
    # DIFFERENTIAL EXPRESSION DATA	
    probes.de = exp.sel
    #	dim(probes.de)
    probes.no.nas = na.omit(dados.anot[match(rownames(probes.de),as.character(dados.anot[,1])),])
    # probes.no.nasAll = na.omit(dados.anot[match(rownames(probes.de),as.character(dados.anot[,1])),])
    #	entrez.remove = rownames(probes.no.nas[duplicated(probes.no.nas[,2]),]) s? repetidos
    probes.no.entrez = probes.no.nas[which(!duplicated(gsub(" *$", "",probes.no.nas[,2]))),]
  }
  
  # ================
  # --- IF AGAIN ---
  # ================
  if("ID" %in% colnames(annotation)){
    annotation = annotation[,c("ID","ENTREZ_GENE_ID","Gene.Symbol")]
    colnames(annotation) = c("probes","entrez","symbol")
    
    # --- REMOVE LINES WITH NULL ---
    sel = annotation[annotation$entrez == "null" | annotation$entrez == "---",] # add other characters <<<<<<<<<<<<<<<<<<<
    if(length(sel[,2]) > 0){
      dados.anot = annotation[-as.numeric(rownames(sel)),] 
    }
    if(length(sel[,2]) == 0){
      dados.anot = annotation
    }
    
    # --- DOUBLE ENTREZ IN SOME LINE ---
    dados.anot[,entrez] = as.character(dados.anot[,entrez])
    dados.anot[,geneSymbol] = as.character(dados.anot[,geneSymbol])
    
    caracter1 = "-"
    caracter2 = "///"
    
    for(i in 1:nrow(dados.anot)){
      tm1 = length(grep(caracter1, dados.anot[i,entrez]))
      tm2 = length(grep(caracter2, dados.anot[i,entrez]))
      
      if(tm1 != 0 || tm2 != 0){
        if(tm1 != 0){
          get1 = strsplit(dados.anot[i,entrez], split=caracter1)[[1]][1]	
          get2 = strsplit(dados.anot[i,geneSymbol], split=caracter1)[[1]][1]	
          dados.anot[i,entrez] = get1
          dados.anot[i,geneSymbol] = get2
        }
        if(tm2 != 0){
          get1 = strsplit(dados.anot[i,entrez], split=caracter2)[[1]][1]
          get2 = strsplit(dados.anot[i,geneSymbol], split=caracter2)[[1]][1]	
          
          dados.anot[i,entrez] = gsub(" *$", "",get1)
          dados.anot[i,geneSymbol] = get2
        }
      }
    }
    expP = expMat[!is.na(match(row.names(expMat),dados.anot$probes)),  ]
    names = dados.anot[match(row.names(expP),dados.anot[,1]),"entrez"]
    row.names(expP) = names
    # OBS: podemos observar aqui a duplicidade de entrezId nas linhas. Selecionei apenas o primeiro entrezId unico, mas futuramente considerar outros metodos.
    expP = expP[which(!duplicated(rownames(expP))),]
    write.table(expP, paste(dirOut,nam.exp,"-MatrixExpressionEntrez.csv",sep=""),sep=",", row.names=T, col.names=T, dec=".", quote=F)
    
    # Matriz diferencialmente expressa
    # exp.sel = matExpProbe[match(rownames(list.DE.all),row.names(matExpProbe)),  ]
    
    exp.sel = expMat[match(list.Exp[[4]],row.names(expMat)),  ]
    # DIFFERENTIAL EXPRESSION DATA	
    probes.de = exp.sel
    probes.no.nas = na.omit(dados.anot[match(rownames(probes.de),as.character(dados.anot[,1])),])
    # probes.no.nasAll = na.omit(dados.anot[match(rownames(probes.de),as.character(dados.anot[,1])),])
    #	entrez.remove = rownames(probes.no.nas[duplicated(probes.no.nas[,2]),]) s? repetidos
    probes.no.entrez = probes.no.nas[which(!duplicated(gsub(" *$", "",probes.no.nas[,2]))),]
  }
  
  # ================
  # --- IF AGAIN ---
  # ================
  if("Gene_ID" %in% colnames(annotation)){
    expP = expMat[!is.na(match(row.names(expMat),as.character(annotation$ID))),  ]
    names = annotation[match(row.names(expP),annotation[,1]),"Gene_ID"]
    row.names(expP) = names
    
    # OBS: podemos observar aqui a duplicidade de entrezId nas linhas. Selecionei apenas o primeiro entrezId unico, mas futuramente considerar outros metodos.
    expP = expP[which(!duplicated(rownames(expP))),]
    write.table(expP, paste(dirOut,nam.exp,"-MatrixExpressionEntrez.csv",sep=""),sep=",", row.names=T, col.names=T, dec=".", quote=F)
    
    exp.sel = expMat[match(list.Exp[[4]],row.names(expMat)),  ]
    # head(exp.sel)
    # DIFFERENTIAL EXPRESSION DATA	
    probes.de = exp.sel
    #	dim(probes.de)
    probes.no.nas = na.omit(annotation[match(rownames(probes.de),as.character(annotation[,1])),])
    # probes.no.nasAll = na.omit(dados.anot[match(rownames(probes.de),as.character(dados.anot[,1])),])
    #	entrez.remove = rownames(probes.no.nas[duplicated(probes.no.nas[,2]),]) s? repetidos
    probes.no.entrez = probes.no.nas[which(!duplicated(gsub(" *$", "",probes.no.nas[,2]))),]   
  }
  
  # =============
  # --- AFTER ---
  # =============
  # MAPPING TO EXPRESSION VALUES 	
  mat.exprs = probes.de[match(probes.no.entrez[,1],rownames(probes.de)),]
  rownames(mat.exprs) = as.character(as.numeric(probes.no.entrez[,2]))
  
  # BIND IN ANNOTATION DATA AND REMOVE DUPLICATES	
  colado = rbind(probes.no.entrez,dados.anot)
  anot.only.entrez = colado[which(!duplicated(colado[,2])),]
  
  # MAPPING TO PVALUE AND FC 
  mat.pvalue = list.Exp[[3]] 
  mat.pvalue.result = mat.pvalue[match(anot.only.entrez[,1],rownames(mat.pvalue)),]
  # change probe to entrez
  rownames(mat.pvalue.result) = anot.only.entrez[,2]
  
  # ADD LIST DATA
  selecionados = anot.only.entrez[which(is.na(match(as.character(anot.only.entrez[,1]),list.Exp[[4]]))==FALSE),]
  list.Exp = list.Exp[-4]
  
  names.list.Exp = c(names(list.Exp),"annotation","selected")
  list.Exp[[2]] = expP
  list.Exp[[3]] =  mat.pvalue.result
  list.Exp[[4]] = anot.only.entrez
  list.Exp[[5]] = selecionados
  names(list.Exp) = names.list.Exp  
  
  # --- create a matrix with data ---
  dadosAllSel = list.Exp$selected
  DE.mat = list.Exp$mat.FC.pVal[match(dadosAllSel$entrez,rownames(list.Exp$mat.FC.pVal)),]
  Exp.mat = exprs(list.Exp[[1]][match(dadosAllSel$entrez,rownames(exprs(list.Exp[[1]]))),])
  dadosAllSel = cbind(dadosAllSel,DE.mat)
  
  write.table(dadosAllSel, paste(dirOut,nam.exp,"-DadosAllSelected.csv",sep=""),sep=",", row.names=F, col.names=T, dec=".", quote=F)
  
  list.Exp
}