# Install required packages if not installed
# install.packages(c("shiny", "shinyjs", "plotly", "ggpointdensity"))

# Load required libraries
library(shiny)
library(shinyjs)
library(plotly)
library(ggpointdensity)


# Define UI
ui <- fluidPage(
  titlePanel("Gene Expression vs ATAC-seq Abundance"),
  sidebarLayout(
    sidebarPanel(
      fileInput("gene_file", "Upload Gene Expression File (TXT format)"),
      fileInput("atac_file", "Upload ATAC-seq Abundance File (TXT format)"),
      radioButtons("plot_type", "Choose Plot Type:",
                   choices = list("Scatter Plot" = "scatter", "Ranked List" = "ranked"),
                   selected = "scatter"),
      actionButton("generate_button", "Generate Figure")

    ),
    mainPanel(
      plotlyOutput("scatter_plot")
    )
  )
)


# atac_raw <- read.table("/Users/peiyu/Documents/GitHub/ATACgraph2/atac.txt")
# gene_raw <- read.table("/Users/peiyu/Documents/GitHub/ATACgraph2/rna.txt")

# Define server
server <- function(input, output) {
  observeEvent(input$generate_button, {
    # Check if both files are uploaded
    if (is.null(input$gene_file) || is.null(input$atac_file)) {
      return(NULL)
    }

    # Read gene expression file
    gene_raw <- read.table(input$gene_file$datapath, header = FALSE)
    gene_raw <- gene_raw[order(gene_raw[,'V2']),] 
     gene_raw$rank <- seq(1, dim(gene_raw)[1])

    # Read ATAC-seq abundance file
    atac_raw <- read.table(input$atac_file$datapath, header = FALSE)

    # Merging gene expression and ATAC-seq data
    data_raw <- merge(gene_raw, atac_raw, by = "V1", all = TRUE)
    colnames(data_raw) <- c("gene", "rna", "rank", "atac")
    data_raw <- na.omit(data_raw)


    # Generating the scatter plot or ranked list based on user choice
    if (input$plot_type == "scatter") {
      output$scatter_plot <- renderPlotly({
          ggplotly(ggplot(aes(y=atac, x=rna), data=data_raw) +  
          # geom_point(alpha = 0.6,colour = "gray")+
          geom_smooth(method = "lm", formula = y~x,se = FALSE,colour="black") +
          theme_classic()+
          geom_pointdensity(alpha = 0.6 ) +
          scale_fill_grey()
          ) %>%
          layout(yaxis = list(title = "ATAC-seq Abundance"),
                 xaxis = list(title = "Gene expression"),
                 showlegend = FALSE,
                annotations = list(x=max(data_raw$rna),
                                    y=max(data_raw$atac),
                                    text = paste0('R',"= ",round(cor(data_raw$rna, data_raw$atac, method = 'spearman'),3),"\n",", p= ", cor.test(data_raw$rna, data_raw$atac, method = 'spearman')$p.value),
                                    showarrow = F
                                    )
                                    )
      })
    } else {
      output$scatter_plot <- renderPlotly({
          ggplotly(ggplot(aes(y=atac, x=rank), data=data_raw) +  
          geom_point(alpha = 0.6,colour = "gray")+
          geom_smooth(method = "lm", formula = y~x,se = FALSE,colour="black") +
          theme_classic()
          )%>%
          layout(yaxis = list(title = "ATAC-seq Abundance"),
                 xaxis = list(title = "Rank of gene expression"),
                 showlegend = FALSE,
                annotations = list(x=max(data_raw$rank),
                                    y=max(data_raw$atac),
                                    text = paste0('R',"= ",round(cor(data_raw$rank, data_raw$atac, method = 'spearman'),3),"\n",", p= ",cor.test(data_raw$rank, data_raw$atac, method = 'spearman')$p.value),
                                    showarrow = F
                                    ))
      })
    }
  })
}

# Run the app
shinyApp(ui, server)



