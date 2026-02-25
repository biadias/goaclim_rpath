library(shiny)
library(Rpath)
library(DT)

setwd("../..")
source("R/EBS_fitting_setup.r")
setwd("apps/ebsfit_1")
run.init <- rsim.fit.run(NA, NA, NA, 
                         scene=scene_fit, run_method="AB", 
                         run_years=hind_years, verbose=T)
fit.init <- rsim.fit.table(scene_fit, run.init)
obs_all     <- strsplit(rsim.fit.list.bio.series(scene_fit)$all,":")
obs_species <- sapply(obs_all, "[[", 1)
obs_sources <- sapply(obs_all, "[[", 2)

ui <- fluidPage(
  # App title ----
  titlePanel("Rpath Fitting Tool"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      selectInput("species","Functional Group",
                  as.character(scene_fit$params$spname)),
      sliderInput(inputId = "pred_v", label = "Pred V",
                  min = -7, max = 7, value = 0, step=0.1),
      sliderInput(inputId = "prey_v", label = "Prey V",
                  min = -7, max = 7, value = 0, step=0.1),  
      plotOutput(outputId = "speciesPlot"),
      ),
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      #plotOutput(outputId = "speciesPlot"),
      div(dataTableOutput(outputId="fitTable"),style="font-size:50%")
      #plotOutput(outputId = "predPlot"),
      #plotOutput(outputId = "preyPlot")
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  #output$speciesControl <- renderUI({
  #  selectInput("speciesChoice", "Functional Group", as.character(scene_fit$params$spname))   
  #})
  run.out <- reactive({
    rsim.fit.run(c(input$pred_v,input$prey_v), 
                 c(input$species,input$species), 
                 c("predvul","preyvul"), 
                 scene=scene_fit, run_method="AB", 
                 run_years=hind_years, verbose=T)  
  })
  
  fit.out <- reactive({
    rsim.fit.table(scene_fit, run.out())
  })
  
  output$speciesPlot <- renderPlot({
    sources <- obs_sources[obs_species==input$species]
    if (length(sources)==0){    
      rsim.plot.fitbio.small(scene_fit, run.out(), input$species, NA)
    } else{
      par(mfrow=c(1,length(sources)))
      for(j in 1:length(sources)){
        rsim.plot.fitbio.small(scene_fit, run.out(), input$species, sources[j])
      }
    }    
  })
  
  output$fitTable <- DT::renderDataTable(fit.out())

  #output$preyPlot <- renderPlot({
  #  rsim.plot(run.base)
  #})
  #
  #output$predPlot <- renderPlot({
  #  rsim.plot(run.base)
  #})
  #bins <- seq(min(x), max(x), length.out = input$bins + 1)
  #
  #hist(x, breaks = bins, col = "#75AADB", border = "white",
  #     xlab = "Waiting time to next eruption (in mins)",
  #     main = "Histogram of waiting times")
  #})
}

shinyApp(ui = ui, server = server)