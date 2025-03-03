library(shiny)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(MASS)

ui <- fluidPage(
    titlePanel("Claret-Bruno Model"),
    sidebarLayout(
        sidebarPanel(
            sliderInput("b0", "b0:", min = 0, max = 10, value = 1),
            sliderInput("kg", "kg:", min = 0, max = 1, value = 0.1),
            sliderInput("nsamples", "Number of samples:", min = 1, max = 1000, value = 100),
            sliderInput("mean_p", "Mean of p:", min = 0, max = 10, value = 1, step = 0.1),
            sliderInput("mean_c", "Mean of c:", min = 0, max = 10, value = 0.1, step = 0.1),
            sliderInput("sd_p", "Standard Deviation of p:", min = 0, max = 10, value = 1, step = 0.1),
            sliderInput("sd_c", "Standard Deviation of c:", min = 0, max = 10, value = 0.1, step = 0.1),
            sliderInput("correlation", "Correlation between p and c:", min = -1, max = 1, value = 0, step = 0.1)
        ),
        mainPanel(
            fixedPanel(
                plotOutput("modelPlot"),
                top = 60, right = 20, width = "80%", height = "80%"
            )
        )
    )
)

server <- function(input, output) {
    output$modelPlot <- renderPlot({
        year <- seq(0, 2, by = 0.05)
        b0 <- input$b0
        kg <- input$kg
        sd_p <- input$sd_p
        sd_c <- input$sd_c
        correlation <- input$correlation
        mean_p <- input$mean_p
        mean_c <- input$mean_c
        
        # Covariance matrix
        cov_matrix <- matrix(
            c(sd_p^2, correlation * sd_p * sd_c, correlation * sd_p * sd_c, sd_c^2), 
            nrow = 2
        )
        
        # Generate 100 samples from the bivariate normal distribution
        nsamples <- input$nsamples
        samples <- mvrnorm(nsamples, mu = c(mean_p, mean_c), Sigma = cov_matrix)
        p <- samples[, 1, drop = FALSE]
        c <- samples[, 2, drop = FALSE]

        ystar_matrix <- b0 * exp(rep(kg, nsamples) %*% t(year) - as.vector(p / c) * (1 - exp(- c %*% t(year)))) # nolint
        data <- data.frame(
            year = rep(year, nsamples), 
            ystar = as.vector(t(ystar_matrix)), 
            sample = rep(1:nsamples, each = length(year))
        )

        ggplot(data, aes(x = year, y = ystar, group = sample)) +
            geom_line(alpha = 0.3) +
            theme_minimal() +
            labs(title = "Claret-Bruno Model", x = "Year", y = "Model Output") +
            ylim(0, 10)
    })
}

shinyApp(ui = ui, server = server)
