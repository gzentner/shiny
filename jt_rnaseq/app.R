library(shiny)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggpubr)

m_vs_c_res <- read.delim("m_vs_c_full_results.tsv", sep = "\t")
mut_vs_het_res <- read.delim("mut_vs_het_full_results.tsv", sep = "\t")

# Define UI ----
ui <- fluidPage(navbarPage(
    title = "ERR RNA-seq",
    
    ## Expression page
    tabPanel("Expression", tabsetPanel(
        ## Volcano plot
        tabPanel("Volcano plot", sidebarLayout(
            sidebarPanel(
                selectInput(
                    "volcanodataset",
                    label = "Experiment",
                    choices = c("Whole animal",
                                "Fat body"),
                    selected = "Whole animal"
                ),
                sliderInput(
                    "volcanopointsize",
                    label = "Point size",
                    min = 0.1,
                    max = 2,
                    value = 0.5
                ),
                sliderInput(
                    "volcanoxlim",
                    label = "X-axis limits",
                    min = -20,
                    max = 20,
                    value = c(-20,20)
                ),
                sliderInput(
                    "volcanoylim",
                    label = "X-axis limits",
                    min = 0,
                    max = 300,
                    value = c(0,300)
                ),
                numericInput(
                    "volcanopadj",
                    label = "Adjusted p-value threshold",
                    value = 0.05
                ),
                numericInput(
                    "volcanofc",
                    label = "log2(fold change) threshold",
                    value = 1
                )
            ),
            mainPanel(plotOutput("volcano"),
                      tableOutput("volcanotable"))
        )),
        ## MA plot
        tabPanel("MA plot", sidebarLayout(
            sidebarPanel(
                selectInput(
                    "madataset",
                    label = "Experiment",
                    choices = c("Whole animal",
                                "Fat body"),
                    selected = "Whole animal"
                ),
                sliderInput(
                    "mapointsize",
                    label = "Point size",
                    min = 0.1,
                    max = 2,
                    value = 0.5
                )
            ),
            mainPanel(plotOutput("maplot"))
        ))
    )),
    
    ## GO page
    tabPanel("Functional enrichment", tabsetPanel(
        ## GO Biological Processes
        tabPanel("GO BP", sidebarLayout(sidebarPanel(
            selectInput(
                "dataset",
                label = "Experiment",
                choices = c("Whole animal",
                            "Fat body"),
                selected = "Whole animal"
            )
        ),
        mainPanel())),
        ## Reactome pathways
        tabPanel("Reactome pathways", sidebarLayout(sidebarPanel(
            selectInput(
                "dataset",
                label = "Experiment",
                choices = c("Whole animal",
                            "Fat body"),
                selected = "Whole animal"
            )
        ),
        mainPanel()))
    ))
))






# Define server logic ----
server <- function(input, output) {
    
    output$volcano <- renderPlot({
        data <- switch(input$volcanodataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(significance = case_when(
                   log2FoldChange <= -input$volcanofc & 
                       padj < input$volcanopadj ~ "downregulated",
                   log2FoldChange >= input$volcanofc & 
                       padj < input$volcanopadj ~ "upregulated",
                   TRUE ~ "not significant")
                   )
        
        ggplot(data, aes(log2FoldChange,-log10(padj))) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            geom_point(aes(color = significance), size = input$volcanopointsize) +
            scale_color_manual(values = c("blue", "gray", "magenta")) +
            geom_vline(xintercept = -input$volcanofc, lty = 2) +
            geom_vline(xintercept = input$volcanofc, lty = 2) +
            geom_hline(yintercept = -log10(input$volcanopadj), lty = 2) +
            xlim(input$volcanoxlim) +
            ylim(input$volcanoylim)
    })
    
    output$volcanotable <- renderTable({
        data <- switch(input$volcanodataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(significance = case_when(
                log2FoldChange <= -input$volcanofc & 
                    padj < input$volcanopadj ~ "downregulated",
                log2FoldChange >= input$volcanofc & 
                    padj < input$volcanopadj ~ "upregulated",
                TRUE ~ "not significant")
            ) %>%
            group_by(significance) %>%
            count
    })
    
    output$maplot <- renderPlot({
        data <- switch(input$madataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        ggmaplot(
            data,
            fdr = 0.05,
            fc = 2,
            size = input$mapointsize,
            genenames = as.vector(data$geneid),
            top = 10,
            font.label = c("bold", 10),
            label.rectangle = T,
            font.legend = "bold",
            font.main = "bold",
            palette = c("magenta", "blue", "gray"),
            ggtheme = theme_classic(),
            select.top.method = "padj"
        )
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)