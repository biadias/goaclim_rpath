library(shiny)
library(Rpath)

# Setup variables
    s.base  <- "../GOA_tiny_base.csv" #system.file("extdata", "GOA_tiny_base.csv", package="Rpath", mustWork=T)
    s.diet  <- "../GOA_tiny_diet.csv" #system.file("extdata", "GOA_tiny_diet.csv", package="Rpath", mustWork=T)
    s.ped   <- "../GOA_tiny_ped.csv" #system.file("extdata", "GOA_tiny_ped.csv",  package="Rpath", mustWork=T)
    s.unbal <- read.rpath.params(s.base, s.diet, s.ped)
    s.bal   <- rpath(s.unbal)
    s.TL <- s.bal$TL; names(s.TL)<-s.bal$Group
    s.B <-  s.bal$Biomass; names(s.B)<-s.bal$Group
    s.Landings <- s.bal$Landings; names(s.Landings)<-s.bal$Group
    s.Discards <- s.bal$Discards; names(s.Discards)<-s.bal$Group
    s.DC <- s.bal$DC
    rownames(s.DC) <- c(s.bal$Group[1:(dim(s.DC)[1]-1)],"Import")
    colnames(s.DC) <- s.bal$Group[1:(dim(s.DC)[2])]


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  # App title ----
  titlePanel("Rpath Shiny Lab 0 - simplified Gulf of Alaska"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      selectInput("species","Functional Group",
         as.character(s.unbal$model$Group[1:s.bal$NUM_LIVING]))
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      plotOutput(outputId = "predPlot"),
      plotOutput(outputId = "preyPlot")
    )
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