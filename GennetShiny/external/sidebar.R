plataformas = c("GPL570: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array",
                "GPL571: [HG-U133A_2] Affymetrix Human Genome U133A 2.0 Array",
                "GPL3535: [Rhesus] Affymetrix Rhesus Macaque Genome Array",
                "GPL96: [HG-U133A] Affymetrix Human Genome U133A Array",
                "GPL6244: [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]",
                "GPL570:	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array",
                "GPL1261:	[Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array",
                "GPL6246:	[MoGene-1_0-st] Affymetrix Mouse Gene 1.0 ST Array [transcript (gene) version]",
                "GPL1355:	[Rat230_2] Affymetrix Rat Genome 230 2.0 Array")

organismAnnot = c("Homo sapiens","Macaca mulatta","Mus musculus","Rattus norvegicus")

normMethod = c("MAS5","RMA","fRMA")

clustMethod = c("average (UPGMA)","complete","median(WPGMC)","centroid(UPGMC)")
distMatrix = c("Euclidian","Pearson","Spearman")

column(4,
       wellPanel(
         fluidRow(
           fileInput('file1', 'Upload Phenodata in CSV File [OPTIONAL]',
                     accept=c('text/csv',
                              'text/comma-separated-values,text/plain',
                              '.csv')),
           # "Note: The input file should have a column with name of archives (SAMPLE_NAME) and ",
           # from", a(href="https://raw.githubusercontent.com/onertipaday/ShinyVolcanoPlot/master/data/example.csv","here"),".",
           tags$hr()
         ),
         fluidRow(
           radioButtons('sep', 'Separator',
                        c(Tab='\t',
                          Comma=','
                        ),
                        selected=',')
         )
       ),
       wellPanel(
         h4("Information about the experiment"),
         fluidRow(
           textInput("expName","Experiment name","")
         ),
         helpText("For example, the name of GEO accession number."),
         tags$hr(),
         fluidRow(
           textInput("overalDesign","Overall design","")
         )
       ),
       wellPanel(
         h4("Normalization Parameters"),
         fluidRow(
           box(width = NULL,
               selectInput("normMeth","Type of Normalization",normMethod,selected=normMethod[1],multiple=FALSE)
           )
         )
       ),
       wellPanel(
         h4("Set Parameters to apply Differential Expression"),
         fluidRow(
           sliderInput("lfcr", "Log2(Fold-Change) in abs:", 
                       0.5,10, value = 1, step=0.1),
           tags$hr(),
           sliderInput("lo", "FDR:", 
                       0.001,0.1, value = 0.05)
         )
       ),
       wellPanel(
         h4("Platform Parameters"),
         fluidRow(
           box(width = NULL,
               selectInput("platSel","Chose Platform",plataformas,selected=plataformas[1],multiple=FALSE)
           )
         ),
         tags$hr(),
         fluidRow(
           box(width = NULL,
               selectInput("orgSel","Chose Organism",organismAnnot,selected=organismAnnot[1],multiple=FALSE)
           )
         )
       ),
       wellPanel(
         h4("Clusterization Parameters"),
         fluidRow(
           box(width = NULL,
               selectInput("clustMet","Hierarchical Method",clustMethod,selected=clustMethod[1],multiple=FALSE)
           )
         ),
         tags$hr(),
         fluidRow(
           box(width = NULL,
               selectInput("dMatrix","Choose the distance matrix",distMatrix,selected=distMatrix[1],multiple=FALSE)
           )
         )
       ),
       # wellPanel(
       #   h4("WGCNA"),
       #   fluidRow(
       #     box(width = NULL,
       #         selectInput("clustMet","Hierarchical Method",clustMethod,selected=clustMethod[1],multiple=FALSE)
       #     )
       #   ),
       #   tags$hr(),
       #   fluidRow(
       #     box(width = NULL,
       #         selectInput("dMatrix","Choose the distance matrix",distMatrix,selected=distMatrix[1],multiple=FALSE)
       #     )
       #   )
       # ),
       wellPanel(
         h4("Execute GeNNet"),
         fluidRow(
           box(width = NULL,
               actionButton("execute","Execute!", icon=icon("gear"))
           )
         )
       )
)
