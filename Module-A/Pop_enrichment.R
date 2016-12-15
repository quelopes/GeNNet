pop_enrichment = function(listData,met,graphDB){
  
  #   matGO = mat.GO
  #   met = method
  #   graphDB = graph
  
  matGO = listData[[7]]
  
  # =============
  # === NEO4J ===
  # =============
  print("--- POP DATABASE WITH NODES BP ---")
  
  # === Checking BP nodes ==
  query = "MATCH (b:BP) RETURN b.GOTerm"
  query = "MATCH (b:GENE) RETURN b.entrez"
  
  resp = cypher(graphDB, query)
  GO.BP = matGO[!duplicated(matGO$GOterm),c(1,2)]
  
  # termsGO = mat.GO[!duplicated(mat.GO[,"GOterm"]),c(1,2)]
  
  # --- ADD GO TERMS ---
  for(n in 1:nrow(GO.BP)){
    getOrCreateNode(graphDB, "BP",GOTerm=GO.BP[n,"GOterm"],Term=GO.BP[n,"term"])
  }
  
  # --- RELATIONSHIPS ---
  matGO[,"cluster"] = paste(matGO[,"cluster"],met,sep="_")
  
  print("Creating relationship BP-->cluster")
  # matGO = matGO[1:5,]
  
  for(m in 1:nrow(matGO)){
  #for(m in 1:10){
    query=paste("MATCH (g:CLUSTER),(p:BP) WHERE g.clustInfo='", matGO[m,"cluster"], "' AND p.GOTerm='",matGO[m,"GOterm"],
                "' CREATE UNIQUE (g)-[r:Was_represented{pValue:'",matGO[m,"Pvalue"],"', NoGOsize:'",matGO[m,"NoGOsize"],
                "', EntrezList:'",matGO[m,"entrez"],"'}]->(p)",sep="") 
    cypher(graphDB, query)
  }
}
