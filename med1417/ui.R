library(shiny)
library(shinythemes)
library(shinyWidgets)
library(colourpicker)
library(tidyverse)
library(viridis)
library(readxl)
library(ggpubr)
library(BiocManager)
library(DT)

options(repos = BiocManager::repositories())
getOption("repos")

ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  title = "Med14/17-AID RNA-seq",
                  
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
                          choices = c("Med14-AID nsRNA",
                                      "Med17-AID nsRNA",
                                      "Med14/17-AID nsRNA",
                                      "Sua7-AID nsRNA",
                                      "Med14-AID total RNA",
                                      "Med17-AID total RNA",
                                      "Med14/17-AID total RNA",
                                      "Sua7-AID total RNA"),
                          selected = "Med14-AID nsRNA"
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
                          choices = c("Med14-AID nsRNA",
                                      "Med17-AID nsRNA",
                                      "Med14/17-AID nsRNA",
                                      "Sua7-AID nsRNA",
                                      "Med14-AID total RNA",
                                      "Med17-AID total RNA",
                                      "Med14/17-AID total RNA",
                                      "Sua7-AID total RNA"),
                          selected = "Med14-AID nsRNA"
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
                          choices = c("Med14-AID nsRNA",
                                      "Med17-AID nsRNA",
                                      "Med14/17-AID nsRNA",
                                      "Sua7-AID nsRNA",
                                      "Med14-AID total RNA",
                                      "Med17-AID total RNA",
                                      "Med14/17-AID total RNA",
                                      "Sua7-AID total RNA"),
                          selected = "Med14-AID nsRNA"
                        )
                      ),
                      mainPanel(DT::dataTableOutput("deg_table"))
                    ))
                  ))
                    
                   
                  ))
                