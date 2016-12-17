lapply(list("shiny","RNeo4j","ggplot2","data.table","shinydashboard","shinythemes","networkD3","igraph","visNetwork","d3heatmap"), function(x) library(x, character.only=T))

options(shiny.usecairo=TRUE)

shinyServer(function(input, output, session){
  
  graph = startGraph("http://localhost:7474/db/data/",username="neo4j",password = "graph")
  way = "/home/rstudio/Data/"
  wayResults = "/home/rstudio/Results/"
  phenoD=dir(way)
  
  example=read.table(paste(way,phenoD[grep(".csv",phenoD)],sep=""),sep=",",head=T)
  
  celArq = dir("/home/rstudio/Data/CEL/")
  
  data = reactive({
    inFile = input$file1
    if(is.null(inFile)) {
      dataframe = example          
    } else {
      dataframe = read.csv(
        inFile$datapath,
        sep=input$sep,
        quote='"',
        stringsAsFactors=FALSE
      )}
  })
  
  # === Platform ===
  plat = reactive({
    if(input$platSel == "GPL570: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array"){
      platnum = "GPL570"
    } else if(input$platSel == "GPL571: [HG-U133A_2] Affymetrix Human Genome U133A 2.0 Array"){
      platnum = "GPL571"
    } else if(input$platSel == "GPL3535: [Rhesus] Affymetrix Rhesus Macaque Genome Array"){
      platnum = "GPL3535"
    }else if(input$platSel == "GPL96: [HG-U133A] Affymetrix Human Genome U133A Array"){
      platnum = "GPL96"
    }else if(input$platSel == "GPL6244: [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]"){
      platnum = "GPL6244"
    }else if(input$platSel == "GPL1355	[Rat230_2] Affymetrix Rat Genome 230 2.0 Array"){
      platnum = "GPL1355"
    }else if(input$platSel == "GPL570	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array"){
      platnum = "GPL570"
    }else if(input$platSel == "GPL1261	[Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array"){
      platnum == "GPL1261"
    }else if(input$platSel == "GPL6246	[MoGene-1_0-st] Affymetrix Mouse Gene 1.0 ST Array [transcript (gene) version]"){
      platnum == "GPL6246"
    }
    
    return(platnum)
  })
  
  org = reactive({
    if(input$orgSel == "Homo sapiens"){
      organismo = "org.Hs.eg.db" 
    }else if(input$orgSel == "Macaca mulatta"){
      organismo = "org.Mmu.eg.db"
    }else if(input$orgSel == "Ratus norvegicus"){
      organismo = "org.Mm.eg.db"
    }else if(input$orgSel == "Mus musculus"){
      organismo = "org.Rn.eg.db"
    }
    return(organismo)
  })
  
  clust = reactive({
    if(input$clustMet == "complete"){
      methodClust="complete"
    } else if(input$clustMet == "average (UPGMA)"){
      methodClust="average"
    }else if(input$clustMet == "median(WPGMC)"){
      methodClust="median"
    }else if(input$clustMet == "centroid(UPGMC)"){
      methodClust="centroid"
    }
    return(methodClust)
  })
  
  distM = reactive({
    if(input$dMatrix == "Euclidian"){
      distMat = "Euclidian"
    }else if(input$dMatrix == "Pearson"){
      distMat = "Pearson"
    }else if(input$dMatrix == "Manhattan"){
      distMat = "Manhattan"
    }
    return(distMat)
  })
  
  normMet = reactive({
    if(input$normMeth == "MAS5"){
      normMeth = "mas5"
    }else if(input$normMeth == "RMA"){
      normMeth = "rma"
    }else if(input$normMeth == "fRMA"){
      normMeth = "frma"
    }
    return(normMeth)
  })
  
  output$tablePheno = renderDataTable({
    dat =  data.frame(data())
  })
  
  output$CELfiles = renderDataTable({
    dat2 = data.frame(celArq)
  })
  
  # === EXECUTE! ===
  v = reactiveValues(data = NULL)
  observeEvent(input$execute, {
    progress = shiny::Progress$new(session,min=1,max=15)
    on.exit(progress$close())
    progress$set(message="Analysis in progress ...") # Wait...and time to have some coffee!
        source("/home/rstudio/Module-A/ModuleProcessingShiny.R")
        v$data = moduleProcessingShiny(input$expName,org(),normMet(),input$overalDesign,plat(),plat(),input$lo,input$lfcr,clust(),distM())
        for (i in 1:15){
          progress$set(value=1)
          Sys.sleep(0.1)
        }
        
  })
  
  output$DBscheme = renderVisNetwork({
    if(is.null(v$data)) return()
    graphStructure = v$data$gdbStructure
    nodes=data.frame(id=unique(as.character(c(graphStructure$This,graphStructure$That))))
    nodes=cbind(nodes,label=nodes$id)
    edges=data.frame(from=graphStructure$This,to=graphStructure$That,label=graphStructure$To,title=graphStructure$To,arrows="to")
    edges=edges[!duplicated(edges$title),]
    visNetwork(nodes,edges)%>% visInteraction(navigationButtons = TRUE)
  })
  
  output$DBmetrics = renderDataTable({
      if(is.null(v$data)) return()
      query2="MATCH (n)-[r]->() RETURN type(r) as relationship ,count(*) as NumRelations"
      graphNumRel = data.frame(cypher(graph,query2))
  })
  
  output$plotHeatMap = renderPlot({
    if(is.null(v$data)) return()
    load(paste(wayResults,"ClutPar.rda",sep=""))
    cols = colorRampPalette(brewer.pal(10, "RdBu"))(256)
    heatmap.3(ClustPar$matExp, Rowv=as.dendrogram(ClustPar$tree),
              scale="none", dendrogram="both", margins=c(6,12),
              Colv=TRUE, ColSideColors=ClustPar$colores,
              symbreaks=FALSE, key=TRUE,
              symkey=FALSE, density.info="none",
              trace="none", #main=main_title, 
              labRow="", cexRow=1, col=rev(cols), ColSideColorsSize=2,
              RowSideColorsSize=2, KeyValueName="z-score")
    legend("topright",legend=c(names(ClustPar$matColor)),
           fill=ClustPar$matColor, border=FALSE, bty="n", y.intersp = 0.9, cex=1)
  #   
  #        # heatmap.2(ClustPar$matExp,trace="none",
  #        #           ColSideColors=ClustPar$colores,
  #        #           density="none", scale="row",
  #        #           cexRow=1.7,cexCol=1.7,
  #        #           labRow="", srtCol=NULL,
  #        #           # offsetCol=0,# offsetRow=0,
  #        #           #margins=c(6,12),
  #        #           Rowv=as.dendrogram(ClustPar$tree)
  #        # )
  #        # legend("topright",legend=c(names(ClustPar$matColor)),
  #        #        fill=ClustPar$matColor, border=FALSE, bty="n", y.intersp = 0.9, cex=1.5)
     })
  
  output$tableFunctionalAnalysis = renderDataTable({
    if(is.null(v$data)) return()
    v$data$functionalAnalises
  })
  
  output$topologies = renderDataTable({
    if(is.null(v$data)) return()
    query4="MATCH (e:EXPERIMENT)-[s:Was_selected]->(g:GENE)-[p:PPI_interaction]-(h:GENE)-[:Was_clusterized]-(c:CLUSTER)-[:Was_represented]-(b:BP) WHERE p.combined_score > 800 RETURN distinct g.symbol,collect(distinct(h.symbol)) AS genes, collect(distinct(b.Term)) AS BP, count(distinct h) AS score order by score DESC Limit 10"
    graphHubs = cypher(graph,query4)
    mat.exp=matrix(0,nrow=nrow(graphHubs),ncol=length(graphHubs))
    colnames(mat.exp) = c("Hub gene", "score", "genes associated", "Biological Process")
    mat.exp[,1] = graphHubs$g.symbol
    mat.exp[,2] = graphHubs$score
    for(i in 1:10){
      mat.exp[i,3] = gsub(",",";",toString(unlist(graphHubs[2][i,1])))
      mat.exp[i,4] =  gsub(",",";",toString(unlist(graphHubs[3][10,1])))
    }
    mat.exp
  })  
  
  output$heatmap <- renderD3heatmap({
    if(is.null(v$data)) return()
    load(paste(wayResults,"ClutPar.rda",sep=""))
    d3heatmap(
      scale(ClustPar$matExp),
      colors = "RdYlBu"#,
      #dendrogram = if (input$cluster) "both" else "none"
    )
  })
  
})
