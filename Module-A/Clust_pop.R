clust_pop = function(listData,usedMethod,graphDB){
  
  # listData = GeNNet
  # graphDB = graph
  # usedMethod = paste(mClust,matDist,sep="-")
  
  nameExp = names(listData[1])
  matClust = listData[[6]]
  
  
  # =============
  # === NEO4J ===
  # =============
  nosUnique = paste("C_",unique(matClust),"_",nameExp,"_",usedMethod,sep="")  
  matClust[1:length(matClust)] = paste("C_",matClust,"_",nameExp,"_",usedMethod,sep="")  
  
  print("--- POP DATABASE WITH NODES CLUSTER ---")
  # === CREATE NODES CLUSTER === 
  
  for(n in 1:length(nosUnique)){
    node = getOrCreateNode(graphDB, "CLUSTER",clustInfo=nosUnique[n])
    query=paste("MATCH (p:CLUSTER),(e:EXPERIMENT) WHERE e.nameExp='", nameExp, "' AND p.clustInfo='",nosUnique[n],
                "' CREATE UNIQUE (p)-[r:Belong]->(e)",sep="") 
    cypher(graphDB, query)
  }
  
  for(i in 1:length(matClust)){
    query=paste("MATCH (g:GENE),(p:CLUSTER) WHERE g.symbol='", names(matClust[i]), "' AND p.clustInfo='",matClust[i],
                "' CREATE UNIQUE (g)-[r:Was_clusterized{Method:'",usedMethod,"'}]->(p)",sep="") 
    cypher(graphDB, query)
  }  
}