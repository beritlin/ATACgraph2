# Install required packages if not installed
# install.packages(c("shiny", "shinyjs", "shinydashboard", "plotly", "ggpointdensity", "shinycssloaders"))

# Load required libraries
library(shiny)
library(shinyjs)
library(shinydashboard)
library(plotly)
library(DT)
library(ggpointdensity)
library(shinycssloaders)
# library(shinythemes)

# Define UI
shinyUI(
  dashboardPage(
    skin = "black",
    # header
    dashboardHeader(title = "ATACgraph2.0"),
    # sidebar
    dashboardSidebar(
      sidebarMenu(
        menuItem("Home", tabName = "home", icon = icon("house")),
        menuItem("Correlation Analyses", tabName = "geneexpress", icon = icon("chart-line")),
        menuItem("Differentially Peaks", tabName = "dep", icon = icon("chart-bar")),
        menuItem("Peak Distribution", tabName = "enrich", icon = icon("chart-column")),
        menuItem("Irreproducible Discovery Rate", tabName = "idr", icon = icon("chart-area"))
      )
    ),
    # main body
    dashboardBody(
      tabItems(
        # home tab content
        tabItem(
          tabName = "home",
          h2("Widgets tab content")
        ),
        # First tab content
        tabItem(
          tabName = "geneexpress",
          h2(" Correlation analyses of ATAC-seqs and RNA-seqs"),
          fluidRow(
            box(fileInput("gene_file", "Upload Gene Expression File (TXT format)"),
              fileInput("atac_file", "Upload ATAC-seq Abundance File (TXT format)"),
              radioButtons("plot_type", "Choose Plot Type:",
                choices = list("Scatter Plot" = "scatter", "Ranked List" = "ranked"),
                selected = "scatter"
              ),
              actionButton("generate_button", "Generate Figure"),
              actionButton("clear_button", "Clear All"),
              height = 500
            ),
            box(shinycssloaders::withSpinner(plotlyOutput("scatter_plot"), color = "grey"), height = 500)
          )
        ),
        # Second tab content
        tabItem(
          tabName = "dep",
          h2("Widgets tab content")
        ),
        # Thirf tab content
        tabItem(
          tabName = "enrich",
          h2("Widgets tab content")
        ),
        # four tab content
        tabItem(
          tabName = "idr",
          h2("Irreproducible Discovery Rate (IDR) Analyses"),
          fluidRow(
            box(fileInput("atac_file1_idr", "Upload Replicate 1 Peaks File (.narrowPeaks)"),
              fileInput("atac_file2_idr", "Upload Replicate 2 Peaks File (.narrowPeaks)"),
              actionButton("generate_button_idr", "Generate Figure"),
              actionButton("clear_button_idr", "Clear All"),
              br(),
              br(),
              h4(shinycssloaders::withSpinner(htmlOutput("idr_para_int"), color = "grey")),
              h4(textOutput("idr_para")),
              h4(htmlOutput("idr_value")),
              height = 500
            ),
            tabBox(
              side = "right", height = 500,
              selected = "Rank",
              tabPanel(
                "Replicate 2",
                plotOutput("idr_plot_r2")
              ),
              tabPanel(
                "Repliacte 1",
                plotOutput("idr_plot_r1")
              ),
              tabPanel(
                "Signal",
                shinycssloaders::withSpinner(plotlyOutput("idr_plot_sig"), color = "grey")
              ),
              tabPanel(
                "Rank",
                shinycssloaders::withSpinner(plotlyOutput("idr_plot_rank"), color = "grey")
              )
            ),
            fluidRow(
              column(
                width = 12,
                box(
                  h4("IDR table"),
                  
                  # Input: Choose dataset ----
                  selectInput("idr_dataset", "Choose a dataset:",
                    choices = c("IDR", "Replicate 1", "Replicate 2"),
                 width = 200),

                  # Button
                  downloadButton("idr_downloadData", "Download"),
              br(),
                            br(),
                  # table
                  DTOutput("idr_table"),
                  width = 12
                )
              )
            )
          )
        )
      )
    )
  )
)
