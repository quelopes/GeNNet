cNodeValues = function(list.genes,graphDB){
  
  # list.genes = GeNNet
  # graphDB=graph
  
  print("--- POP DATABASE WITH DATA AND EXPRESSION ---") 
  
  # === consulta a plataforma utilizada ===
  query=paste("MATCH (e:EXPERIMENT)-[r:Tech_Used]->(p:PLATFORM) WHERE e.nameExp='", names(list.genes[1]), "' RETURN p.platAccess",sep="") 
  plat = cypher(graphDB, query)
  
  matExp = list.genes[[2]]
  matExp = matExp[-which(rownames(matExp) == ""),] # remove null line 
  # matExp=matExp[1:200,]
  namesGenes = row.names(matExp)
  
  col.names = colnames(matExp)
  col.names = paste("['",toString(col.names),"']",sep="")
  col.names = gsub(",","','",col.names)
  col.names = paste(col.names,",",sep="")
  
  vet.all = NULL
  
  for(i in 1:nrow(matExp)){
    vet.all = rbind(vet.all,paste("[",toString(round(matExp[i,],3)),"]",sep=""))
  }
  
  #for(i in 1:50){
  for(i in 1:length(namesGenes)){
    query=paste("MATCH (p:PLATFORM),(g:GENE) WHERE p.platAccess='",plat,"' AND g.entrezId='",namesGenes[i],
                "' CREATE UNIQUE (p)-[r:Has]->(g)",sep="") 
    cypher(graphDB,query)
    query2 = paste("MATCH (e:EXPERIMENT),(g:GENE) WHERE e.nameExp='",names(list.genes[1]),"' AND g.entrezId='",namesGenes[i],
                   "' CREATE UNIQUE (e)-[r:Was_normalized{Method:'","NULL","', threshold:'","NULL","',","times:",col.names,"values:",vet.all[i],"}]->(g)",sep="") 
    cypher(graphDB,query2)
  }
  
  selected = list.genes[[5]][,2]
  for(j in 1:length(selected)){
  #for(j in 1:10){
    query=paste("MATCH (g:GENE),(e:EXPERIMENT) WHERE g.entrezId='", selected[j], "' AND e.nameExp='",names(list.genes[1]),
                "' CREATE UNIQUE (e)-[r:Was_selected{Method:'","FC and pvalue","', threshold:'","1 and 0.05","'}]->(g)",sep="") 
    cypher(graphDB, query)
  }
}
