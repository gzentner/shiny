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
                                   tags$h4("Display options"),
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
                                   tags$h4("Save plot"),
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
                                   tags$h4("Display options"),
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
                                   tags$h4("Save plot"),
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
                                            tags$h4("Display options"),
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
                                                                    "plasma",
                                                                    "viridis"),
                                                        selected = "viridis"),
                                            radioButtons("reactome_grid",
                                                         label = "Show grid",
                                                         choices = c("Yes",
                                                                     "No"),
                                                         selected = "Yes"),
                                            icon = icon("palette")
                                          )),
                                   column(9, 
                                          dropdownButton(
                                            tags$h4("Save plot"),
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
                                             selected = "Count"),
                                 sliderInput("reactome_fontsize",
                                             label = "Font size",
                                             min = 6,
                                             max = 24,
                                             value = 12
                                 ),
                               ),
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
                                     value = 1)
                      ),
                      mainPanel(DT::dataTableOutput("reactome_table"))
                    ))
                  ))
                ))
