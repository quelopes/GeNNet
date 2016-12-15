pre_eSet = function(mat,pheno,annot,wayNames){
  
  # mat = GeNNet
  # pheno = pheno
  # annot = orgAnnot
  # wayNames = nameExp
  
  print("*** ACTIVITY e-SET ***")

  rownames(pheno) = pheno$SAMPLE_NAME
  phenoData = new("AnnotatedDataFrame", data=pheno)
  samples = as.character(pheno$SAMPLE_NAME)
  matSel = mat[,samples]
  
  colnames(matSel) = pheno$SAMPLE_NAME
  eSet = new("ExpressionSet", exprs=matSel, phenoData=phenoData, annotation = annot)
  
  eSet.list = c(eSet)
  
  names(eSet.list) = wayNames
  eSet.list
}
