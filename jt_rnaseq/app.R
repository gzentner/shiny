library(shiny)
library(shinythemes)
library(shinyWidgets)
library(tidyverse)
library(org.Dm.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ggpubr)
library(BiocManager)

options(repos = BiocManager::repositories())
getOption("repos")

m_vs_c_res <- read.delim("m_vs_c_full_results.tsv", sep = "\t")
mut_vs_het_res <- read.delim("mut_vs_het_full_results.tsv", sep = "\t")

# Define UI ----
ui <- fluidPage(theme = shinytheme("cerulean"),
    navbarPage(
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
                numericInput("volcano_padj",
                             label = "Adjusted p-value threshold",
                             value = 0.05),
                numericInput("volcano_fc",
                             label = "log2(fold change) threshold",
                             value = 1),
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
                    value = c(-20, 20)
                ),
                sliderInput(
                    "volcano_ylim",
                    label = "Y-axis limits",
                    min = 0,
                    max = 300,
                    value = c(0, 300)
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
                numericInput("ma_padj",
                             label = "Adjusted p-value threshold",
                             value = 0.05),
                numericInput("ma_fc",
                             label = "log2(fold change) threshold",
                             value = 1),
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
                    value = c(0, 20)
                ),
                sliderInput(
                    "ma_ylim",
                    label = "Y-axis limits",
                    min = -20,
                    max = 20,
                    value = c(-20, 20)
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
        ## Reactome pathways
        tabPanel("Reactome pathways", sidebarLayout(sidebarPanel(
            dropdownButton(
                tags$h3("Save plot"),
                numericInput(inputId = "enrichpathway_width", 
                             label = 'Width', value = 8),
                numericInput(inputId = "enrichpathway_height", 
                             label = 'Height', value = 12),
                textInput(inputId = "enrichpathway_name",
                          label = "File name",
                          value = "reactome.png"),
                downloadButton("reactome_download", "Download plot"),
                icon = icon("save")
            ),
            selectInput(
                "enrichpathway_dataset",
                label = "Experiment",
                choices = c("Whole animal",
                            "Fat body"),
                selected = "Whole animal"),
            radioButtons("enrichpathway_direction",
                         label = "Direction of expression change",
                         choices = c("Upregulated",
                                     "Downregulated"),
                         selected = "Upregulated"),
            numericInput("enrichpathway_padj",
                         label = "Adjusted p-value threshold",
                         value = 0.05),
            numericInput("enrichpathway_fc",
                         label = "log2(fold change) threshold",
                         value = 1),
            numericInput("enrichpathway_ncat",
                         label = "Number of ontologies to display",
                         value = 10),
            selectInput("enrichpathway_plottype",
                        label = "Plot type",
                        choices = c("Dot plot",
                                    "Bar plot"),
                        selected = "Dot plot"),
            selectInput("enrichpathway_x",
                        label = "X-axis value",
                        choices = c("Count",
                                    "GeneRatio"),
                        selected = "Count")),
        mainPanel(plotOutput("reactomepa"))
        ))
    ))
))

# Define server logic ----
server <- function(input, output) {
    
    ## Volcano plot function
    output$volcano <- renderPlot({
        data <- switch(input$volcano_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(
                significance = case_when(
                    log2FoldChange <= -input$volcano_fc &
                        padj < input$volcano_padj ~ "downregulated",
                    log2FoldChange >= input$volcano_fc &
                        padj < input$volcano_padj ~ "upregulated",
                    TRUE ~ "not significant"
                )
            )
        
        ggplot(data, aes(log2FoldChange, -log10(padj))) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            geom_point(aes(color = significance), 
                       size = input$volcano_pointsize) +
            scale_color_manual(values = c("blue", "gray", "magenta")) +
            geom_vline(xintercept = -input$volcano_fc,
                       lty = 2) +
            geom_vline(xintercept = input$volcano_fc, lty = 2) +
            geom_hline(yintercept = -log10(input$volcano_padj),
                       lty = 2) +
            xlim(input$volcano_xlim) +
            ylim(input$volcano_ylim)
    })
    
    ## Volcano plot table function
    output$volcanotable <- renderTable({
        data <- switch(input$volcano_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(
                significance = case_when(
                    log2FoldChange <= -input$volcano_fc &
                        padj < input$volcano_padj ~ "downregulated",
                    log2FoldChange >= input$volcano_fc &
                        padj < input$volcano_padj ~ "upregulated",
                    TRUE ~ "not significant"
                )
            ) %>%
            group_by(significance) %>%
            count
    })
    
    ## MA plot function
    output$maplot <- renderPlot({
        data <- switch(input$ma_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        ggmaplot(
            data,
            fdr = input$ma_padj,
            fc = 2 ^ input$ma_fc,
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
    
    ## Expression table function
    output$deg_table <- DT::renderDataTable({
        data <- switch(input$table_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res) %>%
            dplyr:: select(-change)
        },
        extensions = "Buttons",
        options = list(
            order = list(list(2, "desc"), list(9, "asc")),
            dom = "Bfrtpli",
            buttons = c('copy', 'csv', 'excel', 'print')
        )
)
    
    ## ReactomePA function
    reactome_func <- reactive({
        data <- switch(input$enrichpathway_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(
                significance = case_when(
                    log2FoldChange <= -input$enrichpathway_fc &
                        padj < input$enrichpathway_padj ~ "downregulated",
                    log2FoldChange >= input$enrichpathway_fc &
                        padj < input$enrichpathway_padj ~ "upregulated",
                    TRUE ~ "not significant"
                )
            )
        
        plot_type <- switch(input$enrichpathway_plottype,
                            "Dot plot" = dotplot,
                            "Bar plot" = barplot)
        
        plot_color <- switch(input$enrichpathway_plottype,
                             "Dot plot" = scale_color_viridis_c(),
                             "Bar plot" = scale_fill_viridis_c())
        
        downregulated <- filter(data, significance == "downregulated") %>%
            pull(geneid) %>%
            bitr(., fromType = "SYMBOL", 
                 toType = "ENTREZID", OrgDb = org.Dm.eg.db) %>%
            pull(ENTREZID)
        
        upregulated <- filter(data, significance == "upregulated") %>%
            pull(geneid) %>%
            bitr(., fromType = "SYMBOL", 
                 toType = "ENTREZID", OrgDb = org.Dm.eg.db) %>%
            pull(ENTREZID)
        
        data <- switch(input$enrichpathway_direction,
                       "Upregulated" = upregulated,
                       "Downregulated" = downregulated)
        
        reactomepa <- enrichPathway(data, organism = "fly")
        
        plot_type(reactomepa, 
                  showCategory = input$enrichpathway_ncat,
                  input$enrichpathway_x) + 
            plot_color
    })
        
    output$reactomepa <- renderPlot({reactome_func()})
    
    output$reactome_download <- downloadHandler(
        filename = function() {input$enrichpathway_name}, 
        content = function(file) {ggsave(file, reactome_func(), 
                                         width = input$enrichpathway_width,
                                         height = input$enrichpathway_height)}
    )
}


# Run the app ----
shinyApp(ui = ui, server = server)