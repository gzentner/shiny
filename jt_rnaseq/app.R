library(shiny)
library(shinythemes)
library(shinyWidgets)
library(colourpicker)
library(tidyverse)
library(viridis)
library(org.Dm.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ggpubr)
library(BiocManager)
library(DT)

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
                                fluidRow(
                                    column(3,
                                dropdownButton(
                                    selectInput(
                                        "volcano_theme",
                                        label = "ggplot2 theme",
                                        choices = c("bw",
                                                    "classic",
                                                    "gray",
                                                    "minimal"),
                                        selected = "bw"),
                                    colourInput(
                                        "volcano_up_col",
                                        label = "Upregulated color",
                                        value = "magenta"),
                                    colourInput(
                                        "volcano_down_col",
                                        label = "Downregulated color",
                                        value = "blue"),
                                    colourInput(
                                        "volcano_ns_col",
                                        label = "Unchanged color",
                                        value = "gray"),
                                    radioButtons("volcano_grid",
                                                 label = "Show grid",
                                                 choices = c("Yes",
                                                             "No"),
                                                 selected = "No"),
                                    icon = icon("palette")
                                )),
                                column(9, 
                                dropdownButton(
                                    tags$h3("Save plot"),
                                    numericInput(inputId = "volcano_width", 
                                                 label = 'Width', value = 12),
                                    numericInput(inputId = "volcano_height", 
                                                 label = 'Height', value = 12),
                                    textInput(inputId = "volcano_name",
                                              label = "File name",
                                              value = "volcano.pdf"),
                                    downloadButton("volcano_download", "Download plot"),
                                    icon = icon("save")
                                ))),
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
                                    "volcano_fontsize",
                                    label = "Font size",
                                    min = 6,
                                    max = 24,
                                    value = 12
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
                                fluidRow(
                                    column(3,
                                           dropdownButton(
                                               selectInput(
                                                   "ma_theme",
                                                   label = "ggplot2 theme",
                                                   choices = c("bw",
                                                               "classic",
                                                               "gray",
                                                               "minimal"),
                                                   selected = "bw"),
                                               colourInput(
                                                   "ma_up_col",
                                                   label = "Upregulated color",
                                                   value = "magenta"),
                                               colourInput(
                                                   "ma_down_col",
                                                   label = "Downregulated color",
                                                   value = "blue"),
                                               colourInput(
                                                   "ma_ns_col",
                                                   label = "Unchanged color",
                                                   value = "gray"),
                                               radioButtons("ma_grid",
                                                            label = "Show grid",
                                                            choices = c("Yes",
                                                                        "No"),
                                                            selected = "No"),
                                               icon = icon("palette")
                                           )),
                                    column(9, 
                                           dropdownButton(
                                               tags$h3("Save plot"),
                                               numericInput(inputId = "ma_width", 
                                                            label = 'Width', value = 12),
                                               numericInput(inputId = "ma_height", 
                                                            label = 'Height', value = 12),
                                               textInput(inputId = "ma_name",
                                                         label = "File name",
                                                         value = "maplot.pdf"),
                                               downloadButton("ma_download", "Download plot"),
                                               icon = icon("save")
                                           ))),
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
                                    "ma_fontsize",
                                    label = "Font size",
                                    min = 6,
                                    max = 24,
                                    value = 12
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
                            mainPanel(DT::dataTableOutput("deg_table"))
                        ))
                    )),
                    
                    ## GO page
                    tabPanel("Functional enrichment", tabsetPanel(
                        ## Reactome pathways
                        tabPanel("Reactome pathways", 
                                 sidebarLayout(
                                     sidebarPanel(
                                         fluidRow(
                                             column(3,
                                                    dropdownButton(
                                                        selectInput(
                                                            "reactome_theme",
                                                            label = "ggplot2 theme",
                                                            choices = c("bw",
                                                                        "classic",
                                                                        "gray",
                                                                        "minimal"),
                                                            selected = "bw"),
                                                        selectInput("reactome_pal",
                                                                    label = "Palette",
                                                                    choices = c("cividis",
                                                                                "inferno",
                                                                                "magma",
                                                                                "viridis"),
                                                                    selected = "viridis"),
                                                        radioButtons("reactome_grid",
                                                                     label = "Show grid",
                                                                     choices = c("Yes",
                                                                                 "No"),
                                                                     selected = "No"),
                                                        icon = icon("palette")
                                                    )),
                                             column(9, 
                                                    dropdownButton(
                                                        tags$h3("Save plot"),
                                                        numericInput(inputId = "reactome_width", 
                                                                     label = 'Width', value = 12),
                                                        numericInput(inputId = "reactome_height", 
                                                                     label = 'Height', value = 12),
                                                        textInput(inputId = "reactome_name",
                                                                  label = "File name",
                                                                  value = "reactome.pdf"),
                                                        downloadButton("reactome_download", "Download plot"),
                                                        icon = icon("save")
                                                    ))),
                                         selectInput(
                                             "reactome_dataset",
                                             label = "Experiment",
                                             choices = c("Whole animal",
                                                         "Fat body"),
                                             selected = "Whole animal"),
                                         radioButtons("reactome_direction",
                                                      label = "Direction of expression change",
                                                      choices = c("Upregulated",
                                                                  "Downregulated"),
                                                      selected = "Upregulated"),
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
                                                     selected = "Count")),
                                     mainPanel(plotOutput("reactomepa"))
                                 )),
                        ## Table
                        tabPanel("Table", sidebarLayout(
                            sidebarPanel(
                                selectInput(
                                    "reactome_table_dataset",
                                    label = "Experiment",
                                    choices = c("Whole animal",
                                                "Fat body"),
                                    selected = "Whole animal"
                                ),
                                radioButtons("reactome_table_direction",
                                             label = "Direction of expression change",
                                             choices = c("Upregulated",
                                                         "Downregulated"),
                                             selected = "Upregulated"),
                                numericInput("reactome_table_padj",
                                             label = "Adjusted p-value threshold",
                                             value = 0.05),
                                numericInput("reactome_table_fc",
                                             label = "log2(fold change) threshold",
                                             value = 1),
                            ),
                            mainPanel(DT::dataTableOutput("reactome_table"))
                        ))
                    ))
                ))

# Define server logic ----
server <- function(input, output) {
    
    ## Volcano plot function
    volcano_func <- reactive({
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
        
        theme <- switch(input$volcano_theme,
                        "bw" = theme_bw(),
                        "classic" = theme_classic(),
                        "gray" = theme_gray(),
                        "minimal" = theme_minimal())
        
        grid <- switch(input$volcano_grid,
                       "Yes" = element_line(),
                       "No" = element_blank())
        
        ggplot(data, aes(log2FoldChange, -log10(padj))) +
            theme +
            theme(panel.grid = grid,
                text = element_text(size = input$volcano_fontsize)) +
            geom_point(aes(color = significance), 
                       size = input$volcano_pointsize) +
            scale_color_manual(values = c(input$volcano_down_col,
                                          input$volcano_ns_col,
                                          input$volcano_up_col)) +
            geom_vline(xintercept = -input$volcano_fc,
                       lty = 2) +
            geom_vline(xintercept = input$volcano_fc, lty = 2) +
            geom_hline(yintercept = -log10(input$volcano_padj),
                       lty = 2) +
            xlim(input$volcano_xlim) +
            ylim(input$volcano_ylim)
    })
    
    output$volcano <- renderPlot({volcano_func()})
    
    output$volcano_download <- downloadHandler(
        filename = function() {input$volcano_name}, 
        content = function(file) {ggsave(file, volcano_func(), 
                                         width = input$volcano_width,
                                         height = input$volcano_height)}
    )
    
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
    
    ma_func <- reactive({
        data <- switch(input$ma_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res)
        
        theme <- switch(input$ma_theme,
                        "bw" = theme_bw(),
                        "classic" = theme_classic(),
                        "gray" = theme_gray(),
                        "minimal" = theme_minimal())
        
        grid <- switch(input$ma_grid,
                       "Yes" = grids(linetype = "solid"),
                       "No" = grids(linetype = "blank"))
        
        p <- ggmaplot(
            data,
            fdr = input$ma_padj,
            fc = 2 ^ input$ma_fc,
            size = input$ma_pointsize,
            genenames = as.vector(data$geneid),
            top = 10,
            font.label = input$ma_fontsize,
            label.rectangle = T,
            font.legend = input$ma_fontsize,
            font.main = input$ma_fontsize,
            palette = c(input$ma_up_col, input$ma_down_col, input$ma_ns_col),
            ggtheme = theme,
            select.top.method = "padj",
            xlim = input$ma_xlim,
            ylim = input$ma_ylim,
            font.tickslab = input$ma_fontsize,
            font.x = input$ma_fontsize,
            font.y = input$ma_fontsize
        )
            
            p + grid
    })
    
    output$maplot <- renderPlot({ma_func()})
    
    output$ma_download <- downloadHandler(
        filename = function() {input$ma_name}, 
        content = function(file) {ggsave(file, ma_func(), 
                                         width = input$ma_width,
                                         height = input$ma_height)}
    )
    
    ## Expression table function
    output$deg_table <- DT::renderDataTable({
        data <- switch(input$table_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res) %>%
            dplyr::select(-change)
        
        return(data)},
        
        extensions = "Buttons",
        options = list(
            order = list(list(2, "desc"), list(7, "asc")),
            dom = "Bfrtpli",
            buttons = c('copy', 'csv', 'excel', 'print')
        )
    )
    
    ## ReactomePA function
    reactome_func <- reactive({
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
        
        plot_pal <- reactive({switch(input$reactome_pal,
                           "cividis" = "cividis",
                           "inferno" = "inferno",
                           "magma" = "magma",
                           "viridis" = "viridis")})
        
        plot_color <- reactive({
            if(input$reactome_plottype == "Dot plot")
                return(noquote(paste0("scale_color_", plot_pal(), "_c()")))
            else if(input$reactome_plottype == "Bar plot")
                return(noquote(paste0("scale_fill_", plot_pal(), "_c()")))
        })
            
            
            # switch(input$reactome_plottype,
            #                  "Dot plot" = get(noquote(paste0("scale_color_", plot_pal, "_c()"))),
            #                  "Bar plot" = get(noquote(paste0("scale_fill_", plot_pal, "_c()"))))
        
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
        
        data <- switch(input$reactome_direction,
                       "Upregulated" = upregulated,
                       "Downregulated" = downregulated)
        
        reactomepa <- enrichPathway(data, organism = "fly")
        
        plot_type(reactomepa, 
                  showCategory = input$reactome_ncat,
                  input$reactome_x) + 
            plot_color()

    })
    
    output$reactomepa <- renderPlot({reactome_func()})
    
    output$reactome_download <- downloadHandler(
        filename = function() {input$reactome_name}, 
        content = function(file) {ggsave(file, reactome_func(), 
                                         width = input$reactome_width,
                                         height = input$reactome_height)}
    )
    
    ## Reactome table function
    output$reactome_table <- DT::renderDataTable({
        data <- switch(input$reactome_table_dataset,
                       "Whole animal" = m_vs_c_res,
                       "Fat body" = mut_vs_het_res) %>%
            dplyr::select(-change)
        
        data <- data %>%
            mutate(
                significance = case_when(
                    log2FoldChange <= -input$reactome_table_fc &
                        padj < input$reactome_table_padj ~ "downregulated",
                    log2FoldChange >= input$reactome_table_fc &
                        padj < input$reactome_table_padj ~ "upregulated",
                    TRUE ~ "not significant"
                )
            )
        
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
        
        data <- switch(input$reactome_table_direction,
                       "Upregulated" = upregulated,
                       "Downregulated" = downregulated)
        
        reactomepa <- enrichPathway(data, 
                                    organism = "fly", 
                                    pAdjustMethod = "fdr",
                                    qvalueCutoff = 0.05,
                                    readable = T)
        
        return(as.data.frame(reactomepa))},
        
        extensions = "Buttons",
        options = list(
            order = list(list(2, "desc"), list(7, "asc")),
            dom = "Bfrtpli",
            buttons = c('copy', 'csv', 'excel', 'print')
        )
    )
    
}


# Run the app ----
shinyApp(ui = ui, server = server)