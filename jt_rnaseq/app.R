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
                    "volcano_dataset",
                    label = "Experiment",
                    choices = c("Whole animal",
                                "Fat body"),
                    selected = "Whole animal"
                ),
                numericInput(
                    "volcano_padj",
                    label = "Adjusted p-value threshold",
                    value = 0.05
                ),
                numericInput(
                    "volcano_fc",
                    label = "log2(fold change) threshold",
                    value = 1
                ),
                sliderInput(
                    "volcano_pointsize",
                    label = "Point size",
                    min = 0.1,
                    max = 2,
                    value = 0.5
                ),
                sliderInput(
                    "volcano_xlim",
                    label = "X-axis limits",
                    min = -20,
                    max = 20,
                    value = c(-20,20)
                ),
                sliderInput(
                    "volcano_ylim",
                    label = "X-axis limits",
                    min = 0,
                    max = 300,
                    value = c(0,300)
                )
            ),
            mainPanel(plotOutput("volcano"),
                      tableOutput("volcanotable"))
        )),
        ## MA plot
        tabPanel("MA plot", sidebarLayout(
            sidebarPanel(
                selectInput(
                    "ma_dataset",
                    label = "Experiment",
                    choices = c("Whole animal",
                                "Fat body"),
                    selected = "Whole animal"
                ),
                numericInput(
                    "ma_padj",
                    label = "Adjusted p-value threshold",
                    value = 0.05
                ),
                numericInput(
                    "ma_fc",
                    label = "log2(fold change) threshold",
                    value = 1
                ),
                sliderInput(
                    "ma_pointsize",
                    label = "Point size",
                    min = 0.1,
                    max = 2,
                    value = 0.5
                ),
                sliderInput(
                    "ma_xlim",
                    label = "X-axis limits",
                    min = 0,
                    max = 20,
                    value = c(0,20)
                ),
                sliderInput(
                    "ma_ylim",
                    label = "X-axis limits",
                    min = -20,
                    max = 20,
                    value = c(-20,20)
                )
            ),
            mainPanel(plotOutput("maplot"))
        )),
        ## Table
        tabPanel("Table", sidebarLayout(
            sidebarPanel(
                selectInput(
                    "table_dataset",
                    label = "Experiment",
                    choices = c("Whole animal",
                                "Fat body"),
                    selected = "Whole animal"
                )
            ),
            mainPanel(dataTableOutput("deg_table"))
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
        data <- switch(input$volcano_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(significance = case_when(
                   log2FoldChange <= -input$volcano_fc & 
                       padj < input$volcano_padj ~ "downregulated",
                   log2FoldChange >= input$volcano_fc & 
                       padj < input$volcano_padj ~ "upregulated",
                   TRUE ~ "not significant")
                   )
        
        ggplot(data, aes(log2FoldChange,-log10(padj))) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            geom_point(aes(color = significance), size = input$volcano_pointsize) +
            scale_color_manual(values = c("blue", "gray", "magenta")) +
            geom_vline(xintercept = -input$volcano_fc, lty = 2) +
            geom_vline(xintercept = input$volcano_fc, lty = 2) +
            geom_hline(yintercept = -log10(input$volcano_padj), lty = 2) +
            xlim(input$volcano_xlim) +
            ylim(input$volcano_ylim)
    })
    
    output$volcanotable <- renderTable({
        data <- switch(input$volcano_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(significance = case_when(
                log2FoldChange <= -input$volcano_fc & 
                    padj < input$volcano_padj ~ "downregulated",
                log2FoldChange >= input$volcano_fc & 
                    padj < input$volcano_padj ~ "upregulated",
                TRUE ~ "not significant")
            ) %>%
            group_by(significance) %>%
            count
    })
    
    output$maplot <- renderPlot({
        data <- switch(input$ma_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        ggmaplot(
            data,
            fdr = input$ma_padj,
            fc = 2^input$ma_fc,
            size = input$ma_pointsize,
            genenames = as.vector(data$geneid),
            top = 10,
            font.label = c("bold", 14),
            label.rectangle = T,
            font.legend = c("bold", 14),
            font.main = "bold",
            palette = c("magenta", "blue", "gray"),
            ggtheme = theme_classic(),
            select.top.method = "padj",
            xlim = input$ma_xlim,
            ylim = input$ma_ylim
        )
    })
    
    output$deg_table <- renderDataTable({
        data <- switch(input$table_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res) %>%
            select(-change)
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)