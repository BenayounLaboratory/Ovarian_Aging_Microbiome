library(shiny)
library(dplyr)
library(beeswarm)
library(outliers)

fluidPage(
  navbarPage(
    "Ovarian Index Score Calculator",
    
    # Introduction Tab
    tabPanel("Introduction",
             mainPanel(
               h2("Ovarian Health Index Calculator", align = 'center'),
               p("This calculator helps determine ovarian health index scores based on multiple parameters:"),
               tags$ul(
                 tags$li("AMH Level"),
                 tags$li("FSH Level"),
                 tags$li("INHBA Level"),
                 tags$li("Combined Follicle Count")
               ),
               p("The score is calculated based on reference values from Young Female (YF) and Estropausal (EF) groups. Please read the accompanying paper for more details."),
               p("A CSV template has been included with the reference batch processing panel. Please visit XXXX for the template."),
               p("A CSV template can also be found with the user-defined batch processing panel. Please visit XXX for the template."),
               p("Of note, the user-defined batch processing CSV requires Group_ID column to contain c_YF and c_EF to define the control values.")
             )
    ),
    
    # Data Input Tab
    tabPanel("Calculate Score",
             sidebarLayout(
               sidebarPanel(
                 textInput("sample_id", "Sample ID", ""),
                 textInput("group_id", "Group ID", ""),
                 numericInput("amh", "AMH Level", value = NA),
                 numericInput("fsh", "FSH Level", value = NA),
                 numericInput("inhba", "INHBA Level", value = NA),
                 numericInput("follicle_count", "Combined Follicle Count", value = NA),
                 actionButton("calculate", "Calculate Score"),
                 hr(),
                 actionButton('delete', 'Clear All History'),
                 hr(),
                 downloadButton("download_results", "Download Results")
               ),
               mainPanel(
                 h3("Results"),
                 tableOutput("score_table"),
                 plotOutput("score_plot")
               )
             )
    ),
    # Batch Upload Tab
    tabPanel("Batch Processing Using Reference Values",
             sidebarLayout(
               sidebarPanel(
                 fileInput("file", "Upload CSV File",
                           accept = ".csv"),
                 actionButton("process_batch", "Process Batch"),
                 hr(),
                 actionButton("batch_delete", "Clear All History"),
                 hr(),
                 downloadButton("download_batch", "Download Batch Results")
               ),
               mainPanel(
                 h3("Batch Results"),
                 tableOutput("batch_table"),
                 plotOutput("batch_plot")
               )
             )
    ),
    # User-defined Batch Upload Tab
    tabPanel("Batch Processing Using User-Defined Values",
             sidebarLayout(
               sidebarPanel(
                 fileInput("user_file", "Upload CSV File",
                           accept = ".csv"),
                 actionButton("user_process_batch", "Process Batch"),
                 hr(),
                 actionButton("user_batch_delete", "Clear All History"),
                 hr(),
                 downloadButton("user_download_batch", "Download Batch Results")
               ),
               mainPanel(
                 h3("Batch Results"),
                 tableOutput("user_batch_table"),
                 plotOutput("user_batch_plot")
               )
             )
    )
  )
)
