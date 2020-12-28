library(shiny)

# Load data tables
m_vs_c_res <- read.delim("m_vs_c_full_results.tsv", sep = "\t")
mut_vs_het_res <- read.delim("mut_vs_het_full_results.tsv", sep = "\t")


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
    
    theme <- switch(input$reactome_theme,
                    "bw" = theme_bw(),
                    "classic" = theme_classic(),
                    "gray" = theme_gray(),
                    "minimal" = theme_minimal())
    
    pal <- switch(input$reactome_pal,
                  "cividis" = "cividis",
                  "inferno" = "inferno",
                  "magma" = "magma",
                  "plasma" = "plasma",
                  "viridis" = "viridis")
    
    pal_type <- switch(input$reactome_plottype,
                       "Dot plot" = scale_color_viridis_c(option = pal),
                       "Bar plot" = scale_fill_viridis_c(option = pal))
    
    grid <- switch(input$reactome_grid,
                   "Yes" = element_line(),
                   "No" = element_blank())
    
    p <- plot_type(reactomepa, 
                   showCategory = input$reactome_ncat,
                   input$reactome_x) +
      pal_type +
      theme + 
      theme(panel.grid = grid,
            text = element_text(size = input$reactome_fontsize))
    
    return(p)
    
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
