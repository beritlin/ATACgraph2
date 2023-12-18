library(shiny)
library(shinyjs)
library(shinydashboard)
library(plotly)
library(ggpointdensity)
library(shinycssloaders)
library(DT)

# library(processx)
 options(shiny.maxRequestSize = 30*1024^2)

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
        )%>% toWebGL()
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
                 ))%>% toWebGL()
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

  output$idr_plot_rank <- renderPlotly(NULL)
  output$idr_plot_sig <- renderPlotly(NULL)
  output$idr_para_int <- renderText(NULL)

  observeEvent(input$generate_button_idr, {

  output$idr_para_int <- renderText(paste("<b>","Result","</b>"))


    # Check if both files are uploaded
    if (is.null(input$atac_file1_idr) || is.null(input$atac_file2_idr)) {
      return(NULL)
    }

    # calculating idr
    system(paste0("/work1/home/peiyu/anaconda3/envs/atac/bin/idr --samples ", input$atac_file1_idr$datapath," ",input$atac_file2_idr$datapath, "  --input-file-type narrowPeak --rank p.value --output-file /work1/home/peiyu/2023_atac/tmp/idr.txt --log-output-file /work1/home/peiyu/2023_atac/tmp/tmp.log"))


    # output idr
    atac_idr_log <- read.csv("/work1/home/peiyu/2023_atac/tmp/tmp.log", header = FALSE, sep=" ") 
    output$idr_value <- renderText(paste0("Number of peaks passing IDR cutoff of 0.05: ",atac_idr_log[4,10],"<b>", " ",atac_idr_log[4,11],"</b>"))
    output$idr_para <- renderText(paste0("Mu=", unlist(strsplit(atac_idr_log[2,4], "\\["))[2], ", Sigma=", atac_idr_log[2,5], ", Rho=", atac_idr_log[2,6], ", Mixture params=", unlist(strsplit(atac_idr_log[2,7], "\\]"))))

    # plot idr
    atac_idr_raw <- read.table("/work1/home/peiyu/2023_atac/tmp/idr.txt",header=F)
    colnames(atac_idr_raw)[1:12] <- c("Chr","Start","End","Name","Score","Strand","Signal","pValue","qValue","Summit","LocalIDR","GlobalIDR")
    atac_idr_raw$idr <- ifelse(atac_idr_raw$Score > 540, "TRUE", "FALSE")

    atac_idr_raw <- atac_idr_raw[order(atac_idr_raw$V15),]
    atac_idr_raw$rank1 <- seq(1, nrow(atac_idr_raw), 1)
    atac_idr_raw <- atac_idr_raw[order(atac_idr_raw$V19),]
    atac_idr_raw$rank2 <- seq(1, nrow(atac_idr_raw), 1)

    atac_idr_raw$box1 <- as.factor(as.numeric(cut_number(atac_idr_raw$rank1, 20)))
    atac_idr_raw$box2 <- as.factor(as.numeric(cut_number(atac_idr_raw$rank2, 20)))

    # plot 
    output$idr_plot_rank <- renderPlotly({
      ggplotly(ggplot(aes(x=log(V15),y=log(V19),color=idr),data=atac_idr_raw) +
                   theme_classic() +
                   geom_point(size=1,alpha=.5)+
                  scale_color_manual(values=c('indianred','gray23'))) %>%
           layout(yaxis = list(title = "<b> Replicate 2 signal (log)<b> "),
                  xaxis = list(title = "<b> Replicate 1 signal (log)<b> "),
                  legend = list(title=list(text='<b> IDR >= 0.05 : </b>'),x = 0.2, y = 1.2, orientation = 'h')
        )%>% toWebGL()
      })

    output$idr_plot_sig <- renderPlotly({
        ggplotly(ggplot(aes(x=rank1,y=rank2,color=idr),data=atac_idr_raw) +
                   theme_classic() +
                   geom_point(size=1,alpha=.5)+
                  scale_color_manual(values=c('indianred','gray23'))) %>%
           layout(yaxis = list(title = "<b> Replicate 2 rank<b> "),
                  xaxis = list(title = "<b> Replicate 1 rank<b> "),
                  legend = list(title=list(text='<b> IDR < 0.05 : </b>'),x = 0.2, y = 1.2, orientation = 'h')
        )%>% toWebGL()
      })

      output$idr_plot_r1 <- renderPlot(
        ggplot(data=atac_idr_raw,aes(x=rank1,y=GlobalIDR,fill=box1)) +
                geom_point(colour="gray65",size=1,alpha=.5)+ 
                geom_boxplot(colour="black",outlier.shape = NA)+
                theme_classic()+
                scale_fill_manual(values=rep("NA",20))+
                theme(legend.position='none',
                      axis.title = element_text(face="bold",size=18))+
                labs(x = 'Repliacte 1 peak rank',
                      y = 'IDR (-log10)')
      )

      output$idr_plot_r2 <- renderPlot(
        ggplot(data=atac_idr_raw,aes(x=rank2,y=GlobalIDR,fill=box2)) +
                geom_point(colour="gray65",size=1,alpha=.5)+ 
                geom_boxplot(colour="black",outlier.shape = NA)+
                theme_classic()+
                scale_fill_manual(values=rep("NA",20))+
                theme(legend.position='none',
                      axis.title = element_text(face="bold",size=18))+
                labs(x = 'Repliacte 1 peak rank',
                      y = 'IDR (-log10)')
      )

      output$idr_table <- DT::renderDT({
          datatable(atac_idr_raw[,1:12]) 
      })


  # Add clear button functionality
  observeEvent(input$clear_button_idr, {
    # shinyjs::reset("atac_file1_idr")
    # shinyjs::reset("atac_file2_idr")
    output$idr_plot_r1 <- renderPlot(NULL)
    output$idr_plot_r2 <- renderPlot(NULL)
    output$idr_plot_sig <- renderPlotly(NULL)
    output$idr_plot_rank <- renderPlotly(NULL)
    output$idr_value <- renderText(NULL)
    output$idr_para <- renderText(NULL)
    output$idr_para_int <- renderText(NULL)
    output$idr_table <- renderDT(NULL)
   })


}
)

})
