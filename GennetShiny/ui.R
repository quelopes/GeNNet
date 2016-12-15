lapply(list("shinydashboard","shiny","RNeo4j","ggplot2","data.table","shinythemes","networkD3","igraph","visNetwork","d3heatmap"), function(x) library(x, character.only=T))

tabPanelAbout = source("external/about.R",local=T)$value
tabPanelCite = source("external/cite.R",local=T)$value

headerPanel_2 = function(title, h, windowTitle=title) {    
  tagList(
    tags$head(tags$title(windowTitle)),
    h(title)
  )
}
shinyUI(fluidPage(theme=shinytheme("cosmo"),
                  headerPanel_2(
                    HTML('GeNNet Platform
			<a href="http://www.dexl.lncc.br/" target="_blank"><img id="stats_logo" align="right" alt="gennet logo" src="./img/gennet.svg" width="12%" /></a>'
                    ), h3, "GeNNet"
                  ),
                  fluidRow(
                    source("external/sidebar.R",local=T)$value,
                    source("external/main.R",local=T)$value#,
                  )
))
