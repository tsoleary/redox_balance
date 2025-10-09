# ------------------------------------------------------------------------------
# Shiny app
# TS O'Leary
# ------------------------------------------------------------------------------

# Load Libraries
library(shiny)
library(tidyverse)

# Load data
dat <- readRDS(here::here("data/processed/metabolomics/dat.rds")) |> 
  unnest(data) |> 
  filter(Region != "QC") |> 
  nest()
results <- readRDS(here::here("output/metabolomics/results.rds"))
#DAMs <- readRDS(here::here("output/metabolomics/DAMs.rds"))

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Exploring the Metabolomics Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          
          # Select metabolite
          selectInput("metabolite", 
                      "Metabolite", 
                      selected = "NADH/NAD+ Ratio",
                      choices = dat$Metabolite),
        
        # Select metabolite
        checkboxInput("plot_type",
                      "Plot individual genotypes",
                      value = FALSE)),
        
        
        
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("abundance_plot"),
           uiOutput("stats")
        )
      
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$abundance_plot <- renderPlot({
      
      if(input$plot_type) {
        # Generate plot with individual genotypes
        if(input$metabolite != "GSSG/GSH-NEM Ratio") {
          dat$data[[which(dat$Metabolite == input$metabolite)]] |>
            ggplot(aes(y = NormIntensity,
                       x = Temperature,
                       fill = Temperature)) +
            geom_boxplot() +
            ggbeeswarm::geom_beeswarm(shape = 21,
                                      size = 4,
                                      cex = 4,
                                      dodge.width = 0.75) +
            scale_fill_manual(values = c("grey", "firebrick")) +
            scale_y_continuous(name = paste(input$metabolite, "\n(normalized intensity)")) +
            expand_limits(y = 0) +
            scale_x_discrete(name = element_blank()) +
            facet_wrap(~Locale) +
            theme_minimal(base_size = 16)
        } else {
          dat$data[[which(dat$Metabolite == input$metabolite)]] |>
            ggplot(aes(y = NormIntensity^(-1),
                       x = Temperature,
                       fill = Temperature)) +
            geom_boxplot() +
            ggbeeswarm::geom_beeswarm(shape = 21,
                                      size = 3,
                                      cex = 2,
                                      dodge.width = 0.75) +
            scale_fill_manual(values = c("grey90", "firebrick")) +
            scale_y_continuous(trans = "log2",
                               name = paste("GSH-NEM/GSSG Ratio", "\n(normalized intensity)"),
                               breaks = breaks,
                               minor_breaks = minor_breaks) +
            theme_minimal(base_size = 16) +
            facet_wrap(~Locale)
        }
        

      } else {
        if (input$metabolite != "GSSG/GSH-NEM Ratio") {
          dat$data[[which(dat$Metabolite == input$metabolite)]] |>
            ggplot(aes(y = NormIntensity,
                       x = Region,
                       fill = Temperature)) +
            geom_boxplot() +
            ggbeeswarm::geom_beeswarm(shape = 21,
                                      size = 5,
                                      cex = 2,
                                      dodge.width = 0.75) +
            scale_fill_manual(values = c("grey", "firebrick")) +
            scale_y_continuous(name = paste(input$metabolite, "\n(normalized intensity)")) +
            expand_limits(y = 0) +
            scale_x_discrete(name = element_blank()) +
            theme_minimal(base_size = 16)
        } else {
          dat$data[[which(dat$Metabolite == input$metabolite)]] |>
            ggplot(aes(y = NormIntensity^(-1),
                       x = Region,
                       fill = Temperature)) +
            geom_boxplot() +
            ggbeeswarm::geom_beeswarm(shape = 21,
                                      size = 3,
                                      cex = 2,
                                      dodge.width = 0.75) +
            scale_fill_manual(values = c("grey90", "firebrick")) +
            scale_y_continuous(trans = "log2",
                               name = paste("GSH-NEM/GSSG Ratio", "\n(normalized intensity)"),
                               breaks = breaks,
                               minor_breaks = minor_breaks) +
            theme_minimal(base_size = 16)
        }

      }
      

    })
    
    output$stats <- renderUI({
      table_html <- results |> 
        filter(Metabolite == input$metabolite) |>  
        select(Metabolite, term, estimate, std.error, df, statistic, p.value, p.value.adj) |> 
        mutate( p.value.adj = ifelse(p.value.adj < 0.05, 
                                     sprintf("<span style='color: red;'>%.4f</span>", p.value.adj), 
                                     sprintf("%.4f", p.value.adj))) |> 
        kableExtra::kable("html", escape = FALSE) |>
        kableExtra::kable_styling()
      
      HTML(table_html)  # Convert to HTML for rendering in Shiny
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
