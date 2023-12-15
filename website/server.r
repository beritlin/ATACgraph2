library(shiny)
library(shinyjs)
library(shinydashboard)
library(plotly)
library(ggpointdensity)
library(shinycssloaders)
# library(processx)

shinyServer(function(input, output, session) {
  output$scatter_plot <- renderPlotly(NULL)

#### gene expression
  observeEvent(input$generate_button, {
    # Check if both files are uploaded
    if (is.null(input$gene_file) || is.null(input$atac_file)) {
      return(NULL)
    }

    # Read gene expression file
    gene_raw <- read.table(input$gene_file$datapath, header = FALSE)
    gene_raw <- gene_raw[order(gene_raw[, 'V2']), ]
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
        ggplotly(ggplot(aes(y = atac, x = rna), data = data_raw) +
                   theme_minimal() +
                   geom_pointdensity() +
                   geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "pink2")) %>%
          layout(yaxis = list(title = "ATAC-seq Abundance"),
                 xaxis = list(title = "Gene expression"),
                 showlegend = FALSE,
                 annotations = list(x = max(data_raw$rna),
                                    y = max(data_raw$atac),
                                    text = paste0('R', " = ", round(cor(data_raw$rna, data_raw$atac, method = 'spearman'), 3),
                                                  "\n", ", p = ", cor.test(data_raw$rna, data_raw$atac, method = 'spearman')$p.value),
                                    showarrow = FALSE
                 )
        )
      })
    } else {
      output$scatter_plot <- renderPlotly({
        ggplotly(ggplot(aes(y = atac, x = rank), data = data_raw) +
                   geom_point(alpha = 0.6, colour = "lightgrey") +
                   geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black") +
                   theme_minimal()) %>%
          layout(yaxis = list(title = "ATAC-seq Abundance"),
                 xaxis = list(title = "Rank of gene expression"),
                 showlegend = FALSE,
                 annotations = list(x = max(data_raw$rank),
                                    y = max(data_raw$atac),
                                    text = paste0('R', " = ", round(cor(data_raw$rank, data_raw$atac, method = 'spearman'), 3),
                                                  "\n", ", p = ", cor.test(data_raw$rank, data_raw$atac, method = 'spearman')$p.value),
                                    showarrow = FALSE
                 ))
      })
    }
  })

  # Add clear button functionality
  observeEvent(input$clear_button, {
    shinyjs::reset("gene_file")
    shinyjs::reset("atac_file")
    output$scatter_plot <- renderPlotly(NULL)
  })




#### idr

  output$idr_plot <- renderPlotly(NULL)

  observeEvent(input$generate_button_idr, {

    output$idr_plot <- renderPlotly(NULL)

    # Check if both files are uploaded
    if (is.null(input$atac_file1_idr) || is.null(input$atac_file2_idr)) {
      return(NULL)
    }

    # calculating idr
    # system(paste0("/work1/home/peiyu/anaconda3/envs/atac/bin/idr --samples ", input$atac_file1_idr$datapath," ",input$atac_file2_idr$datapath, "  --input-file-type narrowPeak --rank p.value --output-file /work1/home/peiyu/2023_atac/tmp/idr.txt --log-output-file /work1/home/peiyu/2023_atac/tmp/tmp.log"))

    # atac_idr_log <- read.csv("/Users/peiyu/Documents/GitHub/ATACgraph2/website/tmp.log",header=F,sep=" ")
    # output$idr_log <- renderText(paste0("Mu=", unlist(strsplit(a[2,4], "\\["))[2], ", Sigma=", a[2,5], ", Rho=", a[2,6], ", Mixture params=", unlist(strsplit(a[2,7], "\\]"))))
    # output$idr_value <- renderText(paste0("Number of peaks passing IDR cutoff of 0.05: ",a[4,10], a[4,11]))
    
    system(paste0("/work1/home/peiyu/anaconda3/envs/atac/bin/idr --samples ", input$atac_file1_idr$datapath," ",input$atac_file2_idr$datapath, "  --input-file-type narrowPeak --rank p.value --output-file /work1/home/peiyu/2023_atac/tmp/idr.txt --log-output-file /work1/home/peiyu/2023_atac/tmp/tmp.log"))


    
    atac_idr_log <- read.csv("/work1/home/peiyu/2023_atac/tmp/tmp.log", header = FALSE, sep=" ") 
    output$idr_value <- renderText(paste0("Number of peaks passing IDR cutoff of 0.05: ",atac_idr_log[4,10], atac_idr_log[4,11]))
    output$idr_para <- renderText(paste0("Mu=", unlist(strsplit(atac_idr_log[2,4], "\\["))[2], ", Sigma=", atac_idr_log[2,5], ", Rho=", atac_idr_log[2,6], ", Mixture params=", unlist(strsplit(atac_idr_log[2,7], "\\]"))))

    atac_idr_raw <- read.table("/work1/home/peiyu/2023_atac/tmp/idr.txt",header=F)
    atac_idr_raw$idr <- ifelse(atac_idr_raw$V5 > 540, "TRUE", "FALSE")


    output$idr_plot <- renderPlotly({
        ggplotly(ggplot(aes(x=log(V15),y=log(V19),color=idr),data=atac_idr_raw) +
                   theme_classic() +
                   geom_point()) %>%
           layout(yaxis = list(title = "Replicate2 signal (log)"),
                  xaxis = list(title = "Replicate1 signal (log)")
        )
      })

    # output$table <- renderDT(iris)
    # atac_idr_log <- gene_raw[order(gene_raw[, 'V2']), ]
  #   gene_raw$rank <- seq(1, dim(gene_raw)[1])

  #   # Read ATAC-seq abundance file
  #   atac_raw <- read.table(input$atac_file$datapath, header = FALSE)

  #   # Merging gene expression and ATAC-seq data
  #   data_raw <- merge(gene_raw, atac_raw, by = "V1", all = TRUE)
  #   colnames(data_raw) <- c("gene", "rna", "rank", "atac")
  #   data_raw <- na.omit(data_raw)

  #   # Generating the scatter plot or ranked list based on user choice
  #   if (input$plot_type == "scatter") {
  #     output$scatter_plot <- renderPlotly({
  #       ggplotly(ggplot(aes(y = atac, x = rna), data = data_raw) +
  #                  theme_minimal() +
  #                  geom_pointdensity() +
  #                  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "pink2")) %>%
  #         layout(yaxis = list(title = "ATAC-seq Abundance"),
  #                xaxis = list(title = "Gene expression"),
  #                showlegend = FALSE,
  #                annotations = list(x = max(data_raw$rna),
  #                                   y = max(data_raw$atac),
  #                                   text = paste0('R', " = ", round(cor(data_raw$rna, data_raw$atac, method = 'spearman'), 3),
  #                                                 "\n", ", p = ", cor.test(data_raw$rna, data_raw$atac, method = 'spearman')$p.value),
  #                                   showarrow = FALSE
  #                )
  #       )
  #     })
  #   } else {
  #     output$scatter_plot <- renderPlotly({
  #       ggplotly(ggplot(aes(y = atac, x = rank), data = data_raw) +
  #                  geom_point(alpha = 0.6, colour = "lightgrey") +
  #                  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black") +
  #                  theme_minimal()) %>%
  #         layout(yaxis = list(title = "ATAC-seq Abundance"),
  #                xaxis = list(title = "Rank of gene expression"),
  #                showlegend = FALSE,
  #                annotations = list(x = max(data_raw$rank),
  #                                   y = max(data_raw$atac),
  #                                   text = paste0('R', " = ", round(cor(data_raw$rank, data_raw$atac, method = 'spearman'), 3),
  #                                                 "\n", ", p = ", cor.test(data_raw$rank, data_raw$atac, method = 'spearman')$p.value),
  #                                   showarrow = FALSE
  #                ))
  #     })
  #   }
  # })

  # # Add clear button functionality
  # observeEvent(input$clear_button, {
  #   shinyjs::reset("gene_file")
  #   shinyjs::reset("atac_file")
  #   output$scatter_plot <- renderPlotly(NULL)
   })
}
)


library(ggplot2)
library(GenomicRanges)

library(valr)
a <- read.table("/Users/peiyu/Documents/GitHub/ATACgraph2/website/ATAC-maize-5000-1_peakcall_peaks.narrowPeak")
b <- read.table("/Users/peiyu/Documents/GitHub/ATACgraph2/website/ATAC-maize-5000-2_peakcall_peaks.narrowPeak")

head(a)


gr1 <- GRanges(
    seqnames = a$V1,
    ranges = IRanges(a$V2,a$V3,names=a$V4),
    score1 = a$V7)
gr2 <- GRanges(
    seqnames = b$V1,
    ranges = IRanges(b$V2,b$V3,names=b$V4),
    score2= b$V7)

m <- findOverlaps(gr1, gr2)
peaks <- gr1[queryHits(m)]

# Add the metadata from gr2
mcols(peaks) <- cbind.data.frame(
    mcols(peaks),
    mcols(gr2[subjectHits(m)]))

peaks <- data.frame((mcols(peaks)))

    mu <- 2.6
    sigma <- 1.3
    rho <- 0.8
    p <- 0.7
  idr.out <- est.IDR(peaks, mu, sigma, rho, p, eps=0.001, max.ite=20)
# select observations exceeding IDR threshold=0.01
    IDR.level <- 0.01
    x.selected <- select.IDR(x, idr.out$IDR, IDR.level)

length(x.selected)/nrow(peaks)


b <- read.table("/Users/peiyu/Documents/GitHub/ATACgraph2/website/idr.txt",header=F)
b$idr <- ifelse(b$V5 > 540, "TRUE", "FALSE")
ggplot(aes(x=log(V15),y=log(V19),color=idr),data=b) +
geom_point()

b <- b[order(b$V15),]
b$rank1 <- seq(1, nrow(b), 1)
b <- b[order(b$V19),]
b$rank2 <- seq(1, nrow(b), 1)

ggplot(aes(x=rank1,y=rank2,color=idr),data=b) +
geom_point()

ggplot(aes(x=rank1,y=V12),data=b) +
geom_point()





print(paste0("Mu=", unlist(strsplit(a[2,4], "\\["))[2], ", Sigma=", a[2,5], ", Rho=", a[2,6], ", Mixture params=", unlist(strsplit(a[2,7], "\\]"))))
print(paste0("Number of peaks passing IDR cutoff of 0.05: ",a[4,10], a[4,11]))

ggplotly(ggplot(aes(x=rank2,y=V12),data=b) +
geom_point())

library(processx)

run("ls")
system.time(run("sleep", "10", timeout = 1, error_on_status = FALSE))
system.time(
  run(
    "sh", c("-c", "for i in 1 2 3 4 5; do echo $i; sleep 1; done"),
    timeout = 2, error_on_status = FALSE
  )
)



files1=("SRR6761057_peaks.narrowPeak")

files2=("SRR6761058_peaks.narrowPeak")
print(paste0("/work1/home/peiyu/anaconda3/envs/atac/bin/idr --samples ", files1," ",files2, "  --input-file-type narrowPeak --rank p.value --output-file /work1/home/peiyu/2023_atac/tmp/idr.txt -log-output-file /work1/home/peiyu/2023_atac/tmp/tmp.log"))


