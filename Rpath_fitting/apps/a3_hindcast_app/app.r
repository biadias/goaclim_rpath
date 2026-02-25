library(shiny)
library(Rpath)

# Setup variables
source("R/a3_hindcast_fitting_base_app_setup.r")


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  # App title ----
  titlePanel("Rpath Shiny Lab 0 - simplified Gulf of Alaska"),
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      plotOutput(outputId = "predPlot"),
      plotOutput(outputId = "preyPlot")
    )

)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  #output$speciesControl <- renderUI({
  #  selectInput("speciesChoice", "Functional Group", as.character(unbal$model$Group))   
  #})
  output$preyPlot <- renderPlot({
     par(mar=c(5,8,4,1)+.1)
     barplot( rev(s.DC[,input$species]), las=1, horiz=T,
              xlab="Proportion in diet",ylab="", 
              main=paste(input$species,"diet ( trophic level",
                   format(s.TL[input$species],digits=3),")"))
     })

  output$predPlot <- renderPlot({
     nmort <- s.DC[input$species,] * 
            (s.bal$QB * s.bal$Biomass)[1:s.bal$NUM_LIVING]/s.B[input$species]
     fmort <- c(s.Landings[input$species]/s.B[input$species],
                s.Discards[input$species]/s.B[input$species])
     names(fmort) <- c("Landings","Discards")
     mort <- c(fmort,nmort)
     par(mar=c(5,8,4,1)+.1)
     barplot(rev(mort), las=1, horiz=T,
             ylab="",xlab="Mortality Rate (year-1)",
             main=paste(input$species,"mortality sources"))
     })
    #bins <- seq(min(x), max(x), length.out = input$bins + 1)
    #
    #hist(x, breaks = bins, col = "#75AADB", border = "white",
    #     xlab = "Waiting time to next eruption (in mins)",
    #     main = "Histogram of waiting times")
    #})
}

shinyApp(ui = ui, server = server)