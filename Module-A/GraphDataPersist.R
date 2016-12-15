graphDataPersist = function(listData,graph,wayTo){
  
  # wayTo = dirW
  
  pathneo= "/usr/local/neo4j-community-3.0.6/data/databases/graph.db/"
  
  print("--- PERSISTENT DATABASE ---")
  
  
  # === Create a dirs with Database and statistics ===  
  dir.list = c(paste(wayTo,c("graph.db","GraphDB-metrics"),sep=""))
  
  v = 0
  for(i in 1:length(dir.list)){
    if (file.exists(dir.list[i])){
      v = 0
      #       print(paste("Ambiente ja existe.",i,sep=""))
    } else {
      dir.create(file.path(dir.list[i]))
      v = 1
    }
  }
  
  query="MATCH (a)-[r]-(b) WHERE labels(a) <>[] AND labels(b) <>[] RETURN DISTINCT head(labels(a)) AS This, type(r) as To, head(labels(b)) AS That limit 100"
  graphStructure = cypher(graph,query)
  write.table(graphStructure,file="Results/GraphDB-metrics/graphStructure.csv",sep=",",row.names = F,col.names = T)
  
  query2="MATCH (n) RETURN DISTINCT labels(n) as vertice, count(labels(n)) as NumNodes"
  grapNumNodes = cypher(graph,query2)
  write.table(grapNumNodes,file="Results/GraphDB-metrics/grapNumNodes.csv",sep=",",row.names = F,col.names = T)
  
  query3= "MATCH (n)-[r]->() RETURN type(r) as relationship ,count(*) as NumRelations"
  graphNumRel = cypher(graph,query3)
  write.table(graphNumRel,file="Results/GraphDB-metrics/graphNumRel.csv",sep=",",row.names = F,col.names = T)
  
  # -------------------------
  # --- Query: Hubs genes ---
  # -------------------------
  query4="MATCH (e:EXPERIMENT)-[s:Was_selected]->(g:GENE)-[p:PPI_interaction]-(h:GENE)-[:Was_clusterized]-(c:CLUSTER)-[:Was_represented]-(b:BP) WHERE p.combined_score > 800 RETURN distinct g.symbol,collect(distinct(h.symbol)) AS genes, collect(distinct(b.Term)) AS BP, count(distinct h) AS score order by score DESC Limit 10"
  graphHubs = cypher(graph,query4)
  if(is.null(graphHubs)==FALSE){
    mat.exp=matrix(0,nrow=nrow(graphHubs),ncol=length(graphHubs))
    colnames(mat.exp) = c("Hub gene", "score", "genes associated", "Biological Process")
    mat.exp[,1] = graphHubs$g.symbol
    mat.exp[,2] = graphHubs$score
    for(i in 1:10){
      mat.exp[i,3] = gsub(",",";",toString(unlist(graphHubs[2][i,1])))
      mat.exp[i,4] =  gsub(",",";",toString(unlist(graphHubs[3][10,1])))
    }
    
    write.table(mat.exp,file="Results/GraphDB-metrics/Hubs.csv",sep=",",row.names = F,col.names = T)
  }
  # === PERSIST THE DATABASE ===
  file.copy(from=pathneo, to=wayTo, 
            overwrite = TRUE, recursive = TRUE, 
            copy.mode = TRUE)
  
  # === To GeNNet ===
  
  names.list = c(names(listData),"gdbStructure")
  listData[[8]] = graphStructure
  names(listData) = names.list
  
  listData
}

