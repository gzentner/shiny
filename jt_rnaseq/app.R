library(shiny)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggsci)

m_vs_c_res <- read.delim("m_vs_c_full_results.tsv", sep = "\t")
mut_vs_het_res <- read.delim("mut_vs_het_full_results.tsv", sep = "\t")

# Define UI ----
ui <- fluidPage( navbarPage(title = "ERR RNA-seq",
                            
                            ## Expression page
                            tabPanel("Expression", tabsetPanel(
                                ## Volcano plot
                                tabPanel("Volcano plot", sidebarLayout(
                                    sidebarPanel( selectInput("dataset",
                                                              label = "Experiment",
                                                              choices = c("Whole animal", 
                                                                          "Fat body"), 
                                                              selected = "Whole animal"),
                                                  sliderInput("pointsize",
                                                              label = "Point size",
                                                              min = 0.1, max = 2, value = 0.5)),
                                    mainPanel(plotOutput("volcano"))
                                )),
                                ## MA plot     
                                tabPanel("MA plot", sidebarLayout(
                                    sidebarPanel( selectInput("dataset",
                                                              label = "Experiment",
                                                              choices = c("Whole animal", 
                                                                          "Fat body"), 
                                                              selected = "Whole animal"),
                                                  sliderInput("pointsize",
                                                              label = "Point size",
                                                              min = 0.1, max = 2, value = 0.5)),
                                    mainPanel(plotOutput("maplot"))
                                ))
                            )),
                            
                            ## GO page
                            tabPanel("Functional enrichment", tabsetPanel(
                                ## GO Biological Processes
                                tabPanel("GO BP", sidebarLayout(
                                    sidebarPanel( selectInput("dataset",
                                                              label = "Experiment",
                                                              choices = c("Whole animal", 
                                                                          "Fat body"), 
                                                              selected = "Whole animal")),
                                    mainPanel()
                                )),
                                ## Reactome pathways     
                                tabPanel("Reactome pathways", sidebarLayout(
                                    sidebarPanel( selectInput("dataset",
                                                              label = "Experiment",
                                                              choices = c("Whole animal", 
                                                                          "Fat body"), 
                                                              selected = "Whole animal")),
                                    mainPanel()
                                ))
                            ))
))






# Define server logic ----
server <- function(input, output) {
    
    output$volcano <- renderPlot({
        
        data <- switch(input$dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        ggplot(data, aes(log2FoldChange, -log10(padj))) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            geom_point(aes(color = change), size = input$pointsize) +
            scale_color_manual(values = c("blue", "gray", "magenta")) +
            geom_vline(xintercept = -1, lty = 2) + 
            geom_vline(xintercept = 1, lty = 2) +
            geom_hline(yintercept = 1.3, lty = 2)
    })
    
    output$maplot <- renderPlot({
        
        data <- switch(input$dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        ggmaplot(data, fdr = 0.05, fc = 2, size = input$pointsize, genenames = as.vector(data$geneid), top = 10,
                 font.label = c("bold", 10), label.rectangle = T, font.legend = "bold", font.main = "bold", 
                 palette = c("magenta", "blue", "gray"), ggtheme = theme_classic(), select.top.method = "padj")
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)