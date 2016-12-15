cNodeExperiment = function(nodesInfo,graphDB){
  
  #nodesInfo = infoExp
  #graphDB = graph

# print(nodesInfo)
  
  print("*** ACTIVITY POP DATABASE WITH EXPERIMENT ***") 
    
    getOrCreateNode(graphDB,"EXPERIMENT",nameExp=as.character(nodesInfo$namesExp),overalDesign=as.character(nodesInfo$overalDesign))
    getOrCreateNode(graphDB,"PLATFORM",platAccess=as.character(nodesInfo$platAccess),platName=as.character(nodesInfo$platName))
    query=paste("MATCH (p:PLATFORM),(e:EXPERIMENT) WHERE p.platAccess='",as.character(nodesInfo$platAccess),"' AND e.nameExp='",as.character(nodesInfo$namesExp),
                 "' CREATE UNIQUE (e)-[r:Tech_Used]->(p)",sep="") 
    cypher(graphDB,query)
    
   
    
# === ADD INFO UNIVERSE ===  
#   universe = toString(list.genes[[3]][,"Entrez.Gene"])
#   names(list.genes[1])
#   query = paste("MATCH (e:EXPERIMENT)-[r:Used]-(p:PLATFORM) WHERE e.NameExp='",names(list.genes[1]),"'SET p.universe='",universe,"'",sep="")
#   
#   cypher(graph,query)
#   
#   # === CREATING RELATIONSHIP ===
#   print("Creating relationship EXPERIMENT--> GENE")
#   # m=1
#   for(m in 1:nrow(list.genes[[5]])){
#     query=paste("MATCH (g:GENE),(p:EXPERIMENT) WHERE g.entrez='", list.genes[[5]][m,"model"], "' AND p.NameExp='",names(list.genes[1]),
#                 "' CREATE UNIQUE (p)-[r:Was_selected{Method:'","FC and pvalue","', threshold:'","1 and 0.05","'}]->(g)",sep="") 
#     cypher(graph, query)

}


