library(shiny)
library(shinythemes)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggpubr)
library(msigdbr)
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
        ## GO Biological Processes
        tabPanel("GO Biological Process", sidebarLayout(sidebarPanel(
            selectInput(
                "go_bp_dataset",
                label = "Experiment",
                choices = c("Whole animal",
                            "Fat body"),
                selected = "Whole animal"
            ),
            numericInput("go_bp_padj",
                         label = "Adjusted p-value threshold",
                         value = 0.05),
            numericInput("go_bp_fc",
                         label = "log2(fold change) threshold",
                         value = 1),
            numericInput("go_bp_ncat",
                         label = "Number of ontologies to display",
                         value = 10),
            selectInput("go_bp_plottype",
                        label = "Plot type",
                        choices = c("Dot plot",
                                    "Bar plot"),
                        selected = "Dot plot"),
            selectInput("go_bp_x",
                        label = "X-axis value",
                        choices = c("Count",
                                    "GeneRatio"),
                        selected = "Count")
        ),
        mainPanel(plotOutput("go_bp_up"),
                  plotOutput("go_bp_down"))
        )),
        ## Reactome pathways
        tabPanel("Reactome pathways", sidebarLayout(sidebarPanel(
            selectInput(
                "reactome_dataset",
                label = "Experiment",
                choices = c("Whole animal",
                            "Fat body"),
                selected = "Whole animal"
            ),
            numericInput("reactome_padj",
                         label = "Adjusted p-value threshold",
                         value = 0.05),
            numericInput("reactome_fc",
                         label = "log2(fold change) threshold",
                         value = 1),
            numericInput("reactome_ncat",
                         label = "Number of ontologies to display",
                         value = 10),
            selectInput("reactome_plottype",
                        label = "Plot type",
                        choices = c("Dot plot",
                                    "Bar plot"),
                        selected = "Dot plot"),
            selectInput("reactome_x",
                        label = "X-axis value",
                        choices = c("Count",
                                    "GeneRatio"),
                        selected = "Count")
        ),
        mainPanel(plotOutput("reactome_up"),
                  plotOutput("reactome_down"))
        
        )),
        ## KEGG pathways
        tabPanel("KEGG pathways", sidebarLayout(sidebarPanel(
            selectInput(
                "kegg_dataset",
                label = "Experiment",
                choices = c("Whole animal",
                            "Fat body"),
                selected = "Whole animal"
            ),
            numericInput("kegg_padj",
                         label = "Adjusted p-value threshold",
                         value = 0.05),
            numericInput("kegg_fc",
                         label = "log2(fold change) threshold",
                         value = 1),
            numericInput("kegg_ncat",
                         label = "Number of ontologies to display",
                         value = 10),
            selectInput("kegg_plottype",
                        label = "Plot type",
                        choices = c("Dot plot",
                                    "Bar plot"),
                        selected = "Dot plot"),
            selectInput("kegg_x",
                        label = "X-axis value",
                        choices = c("Count",
                                    "GeneRatio"),
                        selected = "Count")
        ),
        mainPanel(plotOutput("kegg_up"),
                  plotOutput("kegg_down"))
        
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
    output$deg_table <- renderDataTable({
        data <- switch(input$table_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res) %>%
            select(-change)
    })
    
    ## GO BP upregulated function
    output$go_bp_up <- renderPlot({
        data <- switch(input$go_bp_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(
                significance = case_when(
                    log2FoldChange <= -input$go_bp_fc &
                        padj < input$go_bp_padj ~ "downregulated",
                    log2FoldChange >= input$go_bp_fc &
                        padj < input$go_bp_padj ~ "upregulated",
                    TRUE ~ "not significant"
                )
            )
        
        plot_type <- switch(input$go_bp_plottype,
                            "Dot plot" = dotplot,
                            "Bar plot" = enrichplot::barplot)
        
        plot_color <- switch(input$go_bp_plottype,
                             "Dot plot" = scale_color_viridis_c(),
                             "Bar plot" = scale_fill_viridis_c())
        
        msig_df = msigdbr(species = "Drosophila melanogaster", 
                          category = "C5",
                          subcategory = "BP")
        
        t2g = msig_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
        
        upregulated <- filter(data, significance == "upregulated") %>%
            pull(geneid)
        
        upregulated_go_bp <- enricher(upregulated, TERM2GENE = t2g)
        
        plot_type(upregulated_go_bp, 
                showCategory = input$go_bp_ncat,
                x = input$go_bp_x,
                title = "Upregulated") + 
            plot_color
        })
    
    ## GO BP downregulated function
    output$go_bp_down <- renderPlot({
        data <- switch(input$go_bp_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(
                significance = case_when(
                    log2FoldChange <= -input$go_bp_fc &
                        padj < input$go_bp_padj ~ "downregulated",
                    log2FoldChange >= input$go_bp_fc &
                        padj < input$go_bp_padj ~ "upregulated",
                    TRUE ~ "not significant"
                )
            )
        
        plot_type <- switch(input$go_bp_plottype,
                            "Dot plot" = dotplot,
                            "Bar plot" = enrichplot::barplot)
        
        plot_color <- switch(input$go_bp_plottype,
                             "Dot plot" = scale_color_viridis_c(),
                             "Bar plot" = scale_fill_viridis_c())
        
        msig_df = msigdbr(species = "Drosophila melanogaster", 
                          category = "C5",
                          subcategory = "BP")
        
        t2g = msig_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
        
        downregulated <- filter(data, significance == "downregulated") %>%
            pull(geneid)
        
        downregulated_go_bp <- enricher(downregulated, TERM2GENE = t2g)
        
        plot_type(downregulated_go_bp, 
                showCategory = input$go_bp_ncat,
                x = input$go_bp_x,
                title = "Downregulated") + 
            plot_color
    })
    
    ## Reactome upregulated function
    output$reactome_up <- renderPlot({
        data <- switch(input$reactome_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(
                significance = case_when(
                    log2FoldChange <= -input$reactome_fc &
                        padj < input$reactome_padj ~ "downregulated",
                    log2FoldChange >= input$reactome_fc &
                        padj < input$reactome_padj ~ "upregulated",
                    TRUE ~ "not significant"
                )
            )
        
        plot_type <- switch(input$reactome_plottype,
                            "Dot plot" = dotplot,
                            "Bar plot" = enrichplot::barplot)
        
        plot_color <- switch(input$reactome_plottype,
                             "Dot plot" = scale_color_viridis_c(),
                             "Bar plot" = scale_fill_viridis_c())
        
        msig_df = msigdbr(species = "Drosophila melanogaster", 
                          category = "C2",
                          subcategory = "CP:REACTOME")
        
        t2g = msig_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
        
        upregulated <- filter(data, significance == "upregulated") %>%
            pull(geneid)
        
        upregulated_reactome <- enricher(upregulated, TERM2GENE = t2g)
        
        plot_type(upregulated_reactome, 
                  showCategory = input$reactome_ncat,
                  x = input$reactome_x,
                  title = "Upregulated") + 
            plot_color
    })
    
    ## Reactome downregulated function
    output$reactome_down <- renderPlot({
        data <- switch(input$reactome_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(
                significance = case_when(
                    log2FoldChange <= -input$reactome_fc &
                        padj < input$reactome_padj ~ "downregulated",
                    log2FoldChange >= input$reactome_fc &
                        padj < input$reactome_padj ~ "upregulated",
                    TRUE ~ "not significant"
                )
            )
        
        plot_type <- switch(input$reactome_plottype,
                            "Dot plot" = dotplot,
                            "Bar plot" = barplot)
        
        plot_color <- switch(input$reactome_plottype,
                             "Dot plot" = scale_color_viridis_c(),
                             "Bar plot" = scale_fill_viridis_c())
        
        msig_df = msigdbr(species = "Drosophila melanogaster", 
                          category = "C2",
                          subcategory = "CP:REACTOME")
        
        t2g = msig_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
        
        downregulated <- filter(data, significance == "downregulated") %>%
            pull(geneid)
        
        downregulated_reactome <- enricher(downregulated, TERM2GENE = t2g)
        
        plot_type(downregulated_reactome, 
                  showCategory = input$reactome_ncat,
                  x = input$reactome_x,
                  title = "Downregulated") + 
            plot_color
    })
    
    ## KEGG pathways upregulated function
    output$kegg_up <- renderPlot({
        data <- switch(input$kegg_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(
                significance = case_when(
                    log2FoldChange <= -input$kegg_fc &
                        padj < input$kegg_padj ~ "downregulated",
                    log2FoldChange >= input$kegg_fc &
                        padj < input$kegg_padj ~ "upregulated",
                    TRUE ~ "not significant"
                )
            )
        
        plot_type <- switch(input$kegg_plottype,
                            "Dot plot" = dotplot,
                            "Bar plot" = barplot)
        
        plot_color <- switch(input$kegg_plottype,
                             "Dot plot" = scale_color_viridis_c(),
                             "Bar plot" = scale_fill_viridis_c())
        
        msig_df = msigdbr(species = "Drosophila melanogaster", 
                          category = "C2",
                          subcategory = "CP:KEGG")
        
        t2g = msig_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
        
        upregulated <- filter(data, significance == "upregulated") %>%
            pull(geneid)
        
        upregulated_kegg <- enricher(upregulated, TERM2GENE = t2g)
        
        plot_type(upregulated_kegg, 
                  showCategory = input$kegg_ncat,
                  input$kegg_x,
                  title = "Upregulated") + 
            plot_color
    })
    
    ## KEGG pathways downregulated function
    output$kegg_down <- renderPlot({
        data <- switch(input$kegg_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        data <- data %>%
            mutate(
                significance = case_when(
                    log2FoldChange <= -input$kegg_fc &
                        padj < input$kegg_padj ~ "downregulated",
                    log2FoldChange >= input$kegg_fc &
                        padj < input$kegg_padj ~ "upregulated",
                    TRUE ~ "not significant"
                )
            )
        
        plot_type <- switch(input$kegg_plottype,
                            "Dot plot" = dotplot,
                            "Bar plot" = barplot)
        
        plot_color <- switch(input$kegg_plottype,
                             "Dot plot" = scale_color_viridis_c(),
                             "Bar plot" = scale_fill_viridis_c())
        
        msig_df = msigdbr(species = "Drosophila melanogaster", 
                          category = "C2",
                          subcategory = "CP:KEGG")
        
        t2g = msig_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
        
        downregulated <- filter(data, significance == "downregulated") %>%
            pull(geneid)
        
        downregulated_kegg <- enricher(downregulated, TERM2GENE = t2g)
        
        plot_type(downregulated_kegg, 
                  showCategory = input$kegg_ncat,
                  x = input$kegg_x,
                  title = "Downregulated") + 
            plot_color
    })
    
}


# Run the app ----
shinyApp(ui = ui, server = server)