library(shiny)

# Load data tables
results <- read_xlsx("allRNA-with-geneclasses-donczew2020.xlsx")

allgenes <- results %>%
  select(Ensembl_Gene_ID)

med14ns <- results %>% dplyr::select(Ensembl_Gene_ID, log2FC_Med14ns, adjpvalue_Med14ns, Geneclass, RPGstatus, class) %>%
  rename(log2FoldChange = 2, padj = 3, SAGA_TFIID = 4, RPG = 5, CR_TFIID = 6)
med17ns <- results %>% dplyr::select(Ensembl_Gene_ID, log2FC_Med17ns, adjpvalue_Med17ns, Geneclass, RPGstatus, class) %>%
  rename(log2FoldChange = 2, padj = 3, SAGA_TFIID = 4, RPG = 5, CR_TFIID = 6)
med1417ns <- results %>% dplyr::select(Ensembl_Gene_ID, log2FC_Med1417ns, adjpvalue_Med1417ns, Geneclass, RPGstatus, class) %>%
  rename(log2FoldChange = 2, padj = 3, SAGA_TFIID = 4, RPG = 5, CR_TFIID = 6)
sua7ns <- results %>% dplyr::select(Ensembl_Gene_ID, log2FC_Sua7ns, adjpvalue_Sua7ns, Geneclass, RPGstatus, class) %>%
  rename(log2FoldChange = 2, padj = 3, SAGA_TFIID = 4, RPG = 5, CR_TFIID = 6)

med14tot <- results %>% dplyr::select(Ensembl_Gene_ID, log2FC_Med14tot, adjpvalue_Med14tot, Geneclass, RPGstatus, class) %>%
  rename(log2FoldChange = 2, padj = 3, SAGA_TFIID = 4, RPG = 5, CR_TFIID = 6)
med17tot <- results %>% dplyr::select(Ensembl_Gene_ID, log2FC_Med17tot, adjpvalue_Med17tot, Geneclass, RPGstatus, class) %>%
  rename(log2FoldChange = 2, padj = 3, SAGA_TFIID = 4, RPG = 5, CR_TFIID = 6)
med1417tot <- results %>% dplyr::select(Ensembl_Gene_ID, log2FC_Med1417tot, adjpvalue_Med1417tot, Geneclass, RPGstatus, class) %>%
  rename(log2FoldChange = 2, padj = 3, SAGA_TFIID = 4, RPG = 5, CR_TFIID = 6)
sua7tot <- results %>% dplyr::select(Ensembl_Gene_ID, log2FC_Sua7tot, adjpvalue_Sua7tot, Geneclass, RPGstatus, class) %>%
  rename(log2FoldChange = 2, padj = 3, SAGA_TFIID = 4, RPG = 5, CR_TFIID = 6)

# Define server logic ----
server <- function(input, output, session) {
  
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
    
    ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
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
  
  # Gene classifiation function
  geneclass_func <- reactive({
  all_results <- results %>% 
    dplyr::select(Ensembl_Gene_ID, contains("log2FC")) %>%
    pivot_longer(-Ensembl_Gene_ID, values_to = "log2FC", names_to = "Sample") %>%
    mutate(Sample = factor(Sample, levels=c("log2FC_Sua7tot", "log2FC_Sua7ns",
                                            "log2FC_Med14tot", "log2FC_Med14ns",
                                            "log2FC_Med17tot", "log2FC_Med17ns", 
                                            "log2FC_Med1417tot", "log2FC_Med1417ns")))
  
  cr_tfiid <- results %>%
    mutate(coactivator = case_when(
      class == "CR" ~ "Coactivator-redundant",
      class == "TFIID" ~ "TFIID-dependent",
      TRUE ~ "uncategorized"
    )) %>%
    dplyr::select(Ensembl_Gene_ID, log2FC_Sua7ns, log2FC_Med14ns, log2FC_Med17ns, log2FC_Med1417ns, coactivator) %>%
    filter(coactivator != "uncategorized") %>%
    pivot_longer(-c(Ensembl_Gene_ID, coactivator), values_to = "log2FC", names_to = "Sample") %>%
    mutate(Sample = factor(Sample, levels=c("log2FC_Sua7ns", "log2FC_Med14ns", "log2FC_Med17ns", "log2FC_Med1417ns")))
    
  saga_tfiid <- results %>% 
    mutate(saga_tfiid_class = case_when(
      Geneclass == "SAGA-dominated" & RPGstatus == "NA" ~ "SAGA-dominated",
      Geneclass == "TFIID-dominated" & RPGstatus == "NA" ~ "TFIID-dominated",
      TRUE ~ "uncategorized"
    )) %>%
    dplyr::select(Ensembl_Gene_ID, log2FC_Sua7ns, log2FC_Med14ns, log2FC_Med17ns, log2FC_Med1417ns, saga_tfiid_class) %>%
    filter(saga_tfiid_class != "uncategorized") %>%
    pivot_longer(-c(Ensembl_Gene_ID, saga_tfiid_class), values_to = "log2FC", names_to = "Sample") %>%
    mutate(saga_tfiid_class = factor(saga_tfiid_class, levels=c("SAGA-dominated", "TFIID-dominated"))) %>%
    mutate(Sample = factor(Sample, levels=c("log2FC_Sua7ns", "log2FC_Med14ns", "log2FC_Med17ns", "log2FC_Med1417ns")))
  
    data <- switch(input$geneset,
                   "All genes (total + nsRNA)" = all_results,
                   "Coactivator-redundant/TFIID-dependent genes (nsRNA)" = cr_tfiid,
                   "SAGA/TFIID-dominated genes (nsRNA)" = saga_tfiid)
    
    theme <- switch(input$geneclass_theme,
                    "bw" = theme_bw(),
                    "classic" = theme_classic(),
                    "gray" = theme_gray(),
                    "minimal" = theme_minimal())
    
    grid <- switch(input$geneclass_grid,
                   "Yes" = element_line(),
                   "No" = element_blank())
    
    pal <- switch(input$geneclass_pal,
                  "cividis" = scale_fill_viridis_d(option = "cividis"),
                  "inferno" = scale_fill_viridis_d(option = "inferno"),
                  "magma" = scale_fill_viridis_d(option = "magma"),
                  "plasma" = scale_fill_viridis_d(option = "plasma"),
                  "viridis" = scale_fill_viridis_d(option = "viridis"))
    
    ggplot(data, aes(x = Sample, y = log2FC, fill = Sample)) +
    theme + 
    theme(panel.grid = grid,
          text = element_text(size = input$geneclass_fontsize),
          axis.text.x.bottom = element_text(angle = 45, vjust = 0.95, hjust = 1)) + 
    stat_boxplot(inherit.aes = TRUE, geom = 'errorbar', position = position_dodge(width = 0.75), width = 0.4) + 
    geom_boxplot(inherit.aes = TRUE, position = position_dodge(width=0.75), width = 0.5, notch = TRUE, outlier.size = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    pal
  
  })

  output$geneclass_boxplot <- renderPlot({geneclass_func()})
  
  output$geneclass_download <- downloadHandler(
    filename = function() {input$geneclass_name}, 
    content = function(file) {ggsave(file, geneclass_func(), 
                                     width = input$geneclass_width,
                                     height = input$geneclass_height)}
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
    
    user_data <- if(is.null(input$table_file)){allgenes
      } else {
        read.delim(input$table_file$datapath, header = F,
                 sep = ",", quote = "")
        }
    
    user_data <- as.data.frame(user_data)
    
    data2 <- data %>%
      filter(Ensembl_Gene_ID %in% user_data)
    
    return(data2)

},
    
    extensions = "Buttons",
    options = list(
      order = list(list(3, "asc")),
      dom = "Bfrtpli",
      buttons = c('copy', 'csv', 'excel', 'print')
    )
  )
}
