library(shiny)
library(grid)
library(lattice)
library(raster)
library(png)
library(jpeg)
library(latticeExtra)

source("risk_score.R")


anp <- read.csv("data/input/ANPratios.csv")


ui <- fluidPage(
  headerPanel("Predicting the risk of selecting mammalian-adapted influenza viruses"),
  sidebarPanel(
    h5("Choose ANP32A splice variant fractions:"),
    sliderInput('X1', 'fraction of X1 splice variant (fX1)', value= 0, min=0, max = 1, step=0.01),
    uiOutput("X3"),
    tags$hr(style="border-color: black;"),
    h4("Choose properties of the viral strains:"),
    selectInput('startpercK', 'Start fraction of 627K-virus', choices = c(1,20,40,80), selected = 1),
    selectInput('betaK', 'Infection rate of 627K-virus', choices = c(8.8e-7, 2.7e-6, 1e-5), selected = 2.7e-6), 
    selectInput('betaE', 'Infection rate of 627E-virus', choices = c(8.8e-7, 2.7e-6, 1e-5), selected = 2.7e-6),
    selectInput('deltaK', 'Death rate of 627K infected cells', choices = c(2.6, 4, 6.1), selected = 4), 
    selectInput('deltaE', 'Death rate of 627E infected cells', choices = c(2.6, 4, 6.1), selected = 4),
    selectInput('cK', 'Death rate of 627K-virus', choices = c(2.4, 3, 3.6), selected = 3), 
    selectInput('cE', 'Death rate of 627E-virus', choices = c(2.4, 3, 3.6), selected = 3),
    actionButton(inputId = 'submit', label = 'Create heatmap')
  ),
  mainPanel(
    h2("Heatmap of risk scores for selecting mammalian adaptations in the influenza virus polymerase"),
    plotOutput('plot1'),
    verbatimTextOutput('text'),
    plotOutput('plot2')
  )
)

server <- function(input, output) {

  
  output$X3 <- renderUI({
    sliderInput('X3', 'fraction of X3 splice variant (fX3)', value= 0, min=0, max = 1- input$X1, step=0.01)
  })
  
  risk.matrix <- eventReactive(input$submit, {
    
#    data.rm[which(data.rm[ , "betaK"] == input$betaK & data.rm[ , "betaE"] == input$betaE & data.rm[ , "deltaK"] == input$deltaK & data.rm[ , "deltaE"] == input$deltaE
#                  & data.rm[ , "cE"] == input$cE  & data.rm[ , "cK"] == input$cK & data.rm[ , "percK"] == input$startpercK  
#                  ), c("fX1", "fX3","riskscore")]
    
    read.csv(paste0("data/heatmaprawdata/riskmatrix_betaK",input$betaK,"_deltaK",input$deltaK,"_cK",input$cK,"_betaE",input$betaE,"_deltaE",input$deltaE,"_cE",input$cE, "_percK",input$startpercK,".csv"))
       
    })
  
  output$plot1 <- renderPlot({
    threshold <- get.0lines(risk.matrix())
     levelplot(riskscore ~ fX1*fX3, risk.matrix(), at = seq(-1,1,0.02), col.regions = col.fred(100,1.1),
              panel = function(...){
                panel.levelplot(...)
                grid.points( anp[,"X1"]/100, anp[,"X3"]/100 , pch = 1:11, gp = gpar(cex=1.2))
                grid.points(input$X1, input$X3, pch=23, gp = gpar(fill="yellow", cex=1.4))
               # grid.text(c("Ms", "C", "D", "T", "Q", "Go", "Sn", "Sw", "B", "Mp", "Gu" ), anp[,"X1"]/100, anp[,"X3"]/100 )
                panel.lines(threshold$fX1, threshold$fX3, lwd=1, col= "gray")
              }
    ) 
    }
  )
  output$text <- renderPrint({
    rs <- risk_score(input$X1, input$X3, as.numeric(input$startpercK))
    matrix(c(1-input$X1 - input$X3, rs), nrow = 1, ncol = 2, dimnames = list("", c("fraction of X2 splice variant (fX2)", "riskscore")))
  }
  )
  
  output$plot2 <- renderPlot({
    plot(0:1, 0:1, axes=F, ylab="", xlab="", pch=NA)
    title( "Legend", cex = 1.2, adj=0)
    legend("topleft", 
           legend = c("New splice variant composition", "", dimnames(anp)[[1]],"", "selection of mammalian adapted variants", "selection of avian adapted variants", "separation between mammalian and avian adaptation"), 
           pch = c(23, NA, 1:11, NA, 15, 15, NA), 
           col = c( 1, NA, rep(1,11),NA, 2,  4, "gray"),
           pt.bg = c("yellow"),
           pt.cex = c(1.5, rep(1,16)), 
           lty = c(rep(0,16),1),
           bty = "n"
    )
  }
  )
  
}

shinyApp(ui = ui, server = server)