# Let's write a quick Shiny app to play around with the Claret-Bruno model.

library(shiny)
library(tidyverse)
library(cowplot)
library(ggplot2)

ui <- fluidPage(
    titlePanel("Claret-Bruno Model"),
    sidebarLayout(
        sidebarPanel(
            sliderInput("b0", "b0:", min = 0, max = 10, value = 1),
            sliderInput("kg", "kg:", min = 0, max = 1, value = 0.1),
            sliderInput("p", "p:", min = 0, max = 10, value = 1),
            sliderInput("c", "c:", min = 0, max = 1, value = 0.1)
        ),
        mainPanel(
            plotOutput("modelPlot")
        )
    )
)

server <- function(input, output) {
    output$modelPlot <- renderPlot({
        year <- seq(0, 2, by = 0.01)
        b0 <- input$b0
        kg <- input$kg
        p <- input$p
        c <- input$c
        
        ystar <- b0 * exp(kg * year - (p / c) * (1 - exp(-c * year)))
        
        data <- data.frame(year = year, ystar = ystar)
        
        ggplot(data, aes(x = year, y = ystar)) +
            geom_line() +
            theme_minimal() +
            labs(title = "Claret-Bruno Model", x = "Year", y = "Model Output") +
            ylim(0, 10)
    })
}

shinyApp(ui = ui, server = server)
