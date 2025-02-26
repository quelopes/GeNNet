column(8,
       tabsetPanel(
         tabPanelAbout(),
         tabPanel("CEL Archives Data", dataTableOutput("CELfiles")),
         tabPanel("PhenoData matrix", dataTableOutput("tablePheno")),
         #tabPanel("Console", verbatimTextOutput("console"))
         # tabPanel("PCA-result",dataTableOutput("PCAresults")),
         tabPanel("Heatmap",plotOutput("plotHeatMap",width = "850px", height = "1000px")),
         tabPanel("Interactive Heatmap", d3heatmapOutput("heatmap", width = "800px", height = "800px")),
         tabPanel("Functional Analysis",dataTableOutput("tableFunctionalAnalysis")),
         tabPanel("Topologies",dataTableOutput("topologies")),
         tabPanel("Database scheme", visNetworkOutput("DBscheme",width = "100%", height = "400px")),
         tabPanel("Graph-DB metrics",dataTableOutput("DBmetrics")),
         # tabPanel("Database statistics",),
         tabPanelCite()
       )
)
