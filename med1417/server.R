library(shiny)

# Load data tables
results <- read_xlsx("Mediator-codon-optimality.xlsx")

med14ns <- results %>% dplyr::select(Ensembl.Gene.ID, log2FC_Med14ns, adjpvalue_Med14ns, baseMean) %>%
  rename(log2FoldChange = 2, padj = 3)
med17ns <- results %>% dplyr::select(Ensembl.Gene.ID, log2FC_Med17ns, adjpvalue_Med17ns, baseMean) %>%
  rename(log2FoldChange = 2, padj = 3)
med1417ns <- results %>% dplyr::select(Ensembl.Gene.ID, log2FC_Med1417ns, adjpvalue_Med1417ns, baseMean) %>%
  rename(log2FoldChange = 2, padj = 3)
sua7ns <- results %>% dplyr::select(Ensembl.Gene.ID, log2FC_Sua7ns, adjpvalue_Sua7ns, baseMean) %>%
  rename(log2FoldChange = 2, padj = 3)

med14tot <- results %>% dplyr::select(Ensembl.Gene.ID, log2FC_Med14tot, adjpvalue_Med14tot, baseMean) %>%
  rename(log2FoldChange = 2, padj = 3)
med17tot <- results %>% dplyr::select(Ensembl.Gene.ID, log2FC_Med17tot, adjpvalue_Med17tot, baseMean) %>%
  rename(log2FoldChange = 2, padj = 3)
med1417tot <- results %>% dplyr::select(Ensembl.Gene.ID, log2FC_Med1417tot, adjpvalue_Med1417tot, baseMean) %>%
  rename(log2FoldChange = 2, padj = 3)
sua7tot <- results %>% dplyr::select(Ensembl.Gene.ID, log2FC_Sua7tot, adjpvalue_Sua7tot, baseMean) %>%
  rename(log2FoldChange = 2, padj = 3)

# Define server logic ----
server <- function(input, output) {
  
  ## Volcano plot function
  volcano_func <- reactive({
    data <- switch(input$volcano_dataset,
                   "Med14-AID nsRNA" = med14ns,
                   "Med17-AID nsRNA" = med17ns,
                   "Med14/17-AID nsRNA" = med1417ns,
                   "Sua7-AID nsRNA" = sua7ns,
                   "Med14-AID total RNA" = med14tot,
                   "Med17-AID total RNA" = med17tot,
                   "Med14/17-AID total RNA" = med1417tot,
                   "Sua7-AID total RNA" = sua7tot)
    
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
                   "Med14-AID nsRNA" = med14ns,
                   "Med17-AID nsRNA" = med17ns,
                   "Med14/17-AID nsRNA" = med1417ns,
                   "Sua7-AID nsRNA" = sua7ns,
                   "Med14-AID total RNA" = med14tot,
                   "Med17-AID total RNA" = med17tot,
                   "Med14/17-AID total RNA" = med1417tot,
                   "Sua7-AID total RNA" = sua7tot)
    
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
                   "Med14-AID nsRNA" = med14ns,
                   "Med17-AID nsRNA" = med17ns,
                   "Med14/17-AID nsRNA" = med1417ns,
                   "Sua7-AID nsRNA" = sua7ns,
                   "Med14-AID total RNA" = med14tot,
                   "Med17-AID total RNA" = med17tot,
                   "Med14/17-AID total RNA" = med1417tot,
                   "Sua7-AID total RNA" = sua7tot)
    
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
      genenames = as.vector(data$Ensembl.Gene.ID),
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
                   "Med14-AID nsRNA" = med14ns,
                   "Med17-AID nsRNA" = med17ns,
                   "Med14/17-AID nsRNA" = med1417ns,
                   "Sua7-AID nsRNA" = sua7ns,
                   "Med14-AID total RNA" = med14tot,
                   "Med17-AID total RNA" = med17tot,
                   "Med14/17-AID total RNA" = med1417tot,
                   "Sua7-AID total RNA" = sua7tot) %>%
      select_all()
    
    return(data)},
    
    extensions = "Buttons",
    options = list(
      order = list(list(3, "asc")),
      dom = "Bfrtpli",
      buttons = c('copy', 'csv', 'excel', 'print')
    )
  )
}
