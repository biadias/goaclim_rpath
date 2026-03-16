library(shiny)
library(Rpath)

# Setup variables
    s.base <- "../EBS_condensed_563_base.csv"  #system.file("extdata", "EBS_condensed_563_base.csv", package="Rpath", mustWork=T)
    s.diet <- "../EBS_condensed_563_diet.csv"  #system.file("extdata", "EBS_condensed_563_diet.csv", package="Rpath", mustWork=T)
    s.unbal <- read.rpath.params(s.base, s.diet)
    s.bal   <- rpath(s.unbal)
    s.scene <- rsim.scenario(s.bal, s.unbal, years=1990:2089)
    all.gear <- c("Trawl", "Cod.pots", "Longline", "Crab.pots", "State.fisheries") 
    s.unfished <- adjust.fishing(s.scene, "ForcedEffort", all.gear, sim.year = 1990:2089, value=0)
    s.run   <-    rsim.run(s.unfished, method="AB", years=1990:2089)
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
  titlePanel("Rpath Shiny Lab - Bering Sea"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      sliderInput(inputId = "forageF", label = "Forage Fish F",
                  min = 0, max = 2, value = 0, step=0.1),

      sliderInput(inputId = "densedep", label = "Bird foraging density dependence",
                  min = -7, max = 7, value = 0, step=0.1),

      sliderInput(inputId = "satiate", label = "Bird prey satiation",
                  min = -7, max = 7, value = -7, step=0.1),

      sliderInput(inputId = "switching", label = "Bird prey switching",
                  min = 0.3, max = 3, value = 1, step=0.1),

      sliderInput(inputId = "searchrate", label = "Bird Forage fish search rate",
                  min = 0.2, max = .6, value = 0.4, step=0.02)
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      plotOutput(outputId = "runPlot"),
      ##plotOutput(outputId = "predPlot"),
      plotOutput(outputId = "preyPlot")
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  #output$speciesControl <- renderUI({
  #  selectInput("speciesChoice", "Functional Group", as.character(unbal$model$Group))   
  #})

  output$runPlot  <- renderPlot({
    s.fished <- adjust.fishing(s.unfished, "ForcedFRate", "Forage.fish", 
                             sim.year = 1990:2089, value=input$forageF)
    s.fished <- adjust.rsim.predprey(s.fished,"VV","Birds","all",1+exp(input$densedep))
    s.fished <- adjust.rsim.predprey(s.fished,"DD","Birds","all",1+exp(-input$satiate))
    s.fished <- adjust.rsim.predprey(s.fished,"QQ","Birds","Forage.fish",input$searchrate)
    s.fished <- adjust.rsim.predprey(s.fished,"HandleSwitch","Birds","all",input$switching)
    s.run   <-    rsim.run(s.fished, method="AB", years=1990:2089)
    TMP.rsim.plot.ylim(s.run,c("Birds","Forage.fish","Pollock"),ylim=c(0,2.5))
  })

  output$preyPlot <- renderPlot({
    s.fished <- adjust.fishing(s.unfished, "ForcedFRate", "Forage.fish", 
                             sim.year = 1990:2089, value=input$forageF)
    s.fished <- adjust.rsim.predprey(s.fished,"VV","Birds","all",1+exp(input$densedep))
    s.fished <- adjust.rsim.predprey(s.fished,"DD","Birds","all",1+exp(-input$satiate))
    s.fished <- adjust.rsim.predprey(s.fished,"QQ","Birds","Forage.fish",input$searchrate)
    s.fished <- adjust.rsim.predprey(s.fished,"HandleSwitch","Birds","all",input$switching)
    s.run   <-    rsim.run(s.fished, method="AB", years=1990:2089)

    ind <- get.rsim.predprey.link(s.fished,"Birds","all")
    diet <- s.run$annual_Qlink["2089",ind]/sum(s.run$annual_Qlink["2089",ind])
    names(diet) <- as.character(s.fished$params$spname[s.fished$params$PreyFrom[ind]+1])

     par(mar=c(5,8,4,1)+.1)
     barplot(rev(diet), las=1, horiz=T,
              xlab="Proportion in diet",ylab="", 
              main="Birds Diet")
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