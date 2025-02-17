library(shiny)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(readxl)
library(DT)
# 定义 UI
ui <- fluidPage(
  titlePanel("血浆蛋白污染校正与差异表达分析"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("data_source", "选择数据来源",
                   choices = list("Load experimental data" = "experimental", 
                                  "Load example data" = "example"),
                   selected = "experimental"),
      conditionalPanel(
        condition = "input.data_source == 'experimental'",
        fileInput("data_file", "上传数据文件 (CSV)", accept = ".csv"),
        fileInput("group_file", "上传分组信息文件 (CSV)", accept = ".csv")
      ),
      selectInput(
        "type",
        "校正类型",
        choices = c("all", "erythrocyte", "platelet", "coagulation")
      ),
      selectInput(
        "group1",
        "分组1",
        choices = NULL
      ),
      selectInput(
        "group2",
        "分组2",
        choices = NULL
      ),
      actionButton("run_check", "运行数据检查"),
      actionButton("run_correct", "运行数据校正"),
      actionButton("run_de", "运行差异表达分析")
    ),
    mainPanel(
      tabsetPanel(
tabPanel("原始数据", 
                 conditionalPanel(
                   condition = "output.data_loaded == false",
                   h4("You did not upload your data. Please upload the expression data, or load the example data to check first.")
                 ),
                 conditionalPanel(
                   condition = "output.data_loaded == true",
                   h3("数据文件"),
                   DTOutput("data_table"), 
                   h3("分组信息文件"),
                   DTOutput("group_table")
                 )),  
        tabPanel("数据检查结果", 
                 h3("缺失基因"),
                 verbatimTextOutput("missing_genes"),
                 h3("QC"),
                 tabsetPanel(
                   tabPanel("PCA",
                            plotOutput("pca_pre_plot")),
                   tabPanel("heatmap", 
                            plotOutput("heatmap_pre_plot")),
                   tabPanel("boxplot", 
                            plotOutput("boxplot_pre_plot"))
                 ),
                 # 二级选项卡：污染情况
                 h3("污染marker表达情况"),
                 tabsetPanel(
                   tabPanel("红系污染",
                            DTOutput("data_marker_erythrocyte"),
                            plotOutput("contamination_erythrocyte_plot")),
                   tabPanel("凝血污染", 
                            DTOutput("data_marker_coagulation"),
                            plotOutput("contamination_coagulation_plot")),
                   tabPanel("血小板污染", 
                            DTOutput("data_marker_platelet"),
                            plotOutput("contamination_platelet_plot"))
                 ),
                 h3("样本污染情况"),
                 tabsetPanel(
                   tabPanel("CV",
                            plotOutput("cv_pre_plot")),
                   tabPanel("红系污染",
                            plotOutput("erythrocyte_marker_pre_plot")),
                   tabPanel("凝血污染", 
                            plotOutput("coagulation_marker_pre_plot")),
                   tabPanel("血小板污染", 
                            plotOutput("platelet_marker_pre_plot"))
                 )
        ),
        tabPanel("校正结果", 
                 h3("QC"),
                 tabsetPanel(
                   tabPanel("PCA",
                            plotOutput("pca_post_plot")),
                   tabPanel("heatmap", 
                            plotOutput("heatmap_post_plot")),
                   tabPanel("boxplot", 
                            plotOutput("boxplot_post_plot"))
                 ),
                 h3("污染检查"),
                 tabsetPanel(
                   tabPanel("CV",
                            plotOutput("cv_post_plot")),
                   tabPanel("红系污染",
                            plotOutput("erythrocyte_marker_post_plot")),
                   tabPanel("凝血污染", 
                            plotOutput("coagulation_marker_post_plot")),
                   tabPanel("血小板污染", 
                            plotOutput("platelet_marker_post_plot"))
                   ),
                 h3("矫正后矩阵"),
                 downloadButton("download_data", "下载校正数据"),
                 DTOutput("data_correct_table")
        ),
        tabPanel("差异表达分析", 
                 tabsetPanel(
                             tabPanel("原始数据table",
                                      downloadButton("download_de_pre", "下载差异表达结果"),
                                      DTOutput("result_de_pre_table")),
                             tabPanel("矫正table",
                                      downloadButton("download_de_post", "下载差异表达结果"),
                                      DTOutput("result_de_post_table")
                                      )),
                 tabsetPanel(
                             tabPanel("原始数据",
                                      plotOutput("volc_de_pre")),
                             tabPanel("矫正",
                                      plotOutput("volc_de_post"))
                             )
                 )
        )
      )
    )
  )


# 定义 Server 逻辑 ----
server <- function(input, output, session) {
  ## 加载示例数据 ----
  example_data <- reactive({
    df <- read.csv("./tests/raw_data_aggr.csv")  # 示例数据路径
    rownames(df) <- df[, 1]
    df[, -1, drop = FALSE]
  })
  
  example_group <- reactive({
    df <- read.csv("./tests/group.csv")  # 示例分组路径
    rownames(df) <- df[, 1]
    df
  })
  ## 读取数据 ----
  data <- reactive({
    if (input$data_source == "example") {
      example_data()
    } else {
      req(input$data_file)
      df <- read.csv(input$data_file$datapath)
      rownames(df) <- df[, 1]
      df[, -1, drop = FALSE]
    }
  })
  
  data_group <- reactive({
    if (input$data_source == "example") {
      example_group()
    } else {
      req(input$group_file)
      df <- read.csv(input$group_file$datapath)
      rownames(df) <- df[, 1]
      df
    }
  })
  # 数据加载状态判断
  output$data_loaded <- reactive({
    !is.null(data()) && !is.null(data_group())
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
  
  ## 展示原始数据 ----
  output$data_table <- renderDT({
    req(data())
    datatable(data(), options = list(pageLength = 10))
  })
  
  output$group_table <- renderDT({
    req(data_group())
    datatable(data_group(), options = list(pageLength = 10))
  })
  
  ## 更新分组选择 ----
  observe({
    req(data_group())
    if (!"group" %in% colnames(data_group())) {
      showNotification("分组文件缺少'group'列!", type = "error")
      return()
    }
    updateSelectInput(session, "group1", choices = unique(data_group()$group))
    updateSelectInput(session, "group2", choices = unique(data_group()$group))
  })
  
  ## 数据检查 ----
  result_check <- eventReactive(input$run_check, {
    data_check(data(), data_group())
  })
  ## 数据校正 ----
  ## 开始数据矫正 ----
  result_correct <- eventReactive(input$run_correct, {
    req(result_check())
    data_correct(data = result_check(), type = input$type)
  })
  ## 显示矫正后矩阵 ----
  output$data_correct_table <- renderDT({
    req(result_correct())
    datatable(result_correct()$correct_data, options = list(pageLength = 10))  # 每页显示 10 行
  })
  ## QC(需修改) ----
  ### pre ----
  output$pca_pre_plot <- renderPlot({
    source("./R/modules/QC_PCA.R")
    req(result_check())
    QC_PCA(data = result_check()$rawdata,
           data_group = result_check()$group)
  })
  output$heatmap_pre_plot <- renderPlot({
    source("./R/modules/QC_heatmap.R")
    req(result_check())
    QC_heatmap(data = result_check()$rawdata,
           data_group = result_check()$group)
  })
  output$boxplot_pre_plot <- renderPlot({
    source("./R/modules/QC_boxplot.R")
    req(result_check())
    QC_boxplot(data = result_check()$rawdata,
           data_group = result_check()$group)
  })
  ### post ----
  output$pca_post_plot <- renderPlot({
    source("./R/modules/QC_PCA.R")
    req(result_correct())
    QC_PCA(data = result_correct()$correct_data,
           data_group = result_correct()$group)
  })
  output$heatmap_post_plot <- renderPlot({
    source("./R/modules/QC_heatmap.R")
    req(result_correct())
    QC_heatmap(data = result_correct()$correct_data,
               data_group = result_correct()$group)
  })
  output$boxplot_post_plot <- renderPlot({
    source("./R/modules/QC_boxplot.R")
    req(result_correct())
    QC_boxplot(data = result_correct()$correct_data,
               data_group = result_correct()$group)
  })
  ## 显示缺失基因 ----
  output$missing_genes <- renderPrint({
    req(result_check())
    cat("红系marker缺失基因:\n")
    print(result_check()$missing_genes$missing_erythrocyte_genes)
    cat("凝血marker缺失基因:\n")
    print(result_check()$missing_genes$missing_coagulation_genes)
    cat("血小板marker缺失基因:\n")
    print(result_check()$missing_genes$missing_platelet_genes)
  })
  
  ## 显示污染水平可视化 ----
  ### CV ----
  output$cv_pre_plot <- renderPlot({
    source("./R/modules/get_cv.R")
    req(result_check())
    data <- result_check()
    result_cv_coa <- get_cv(raw_data = data$rawdata,protein = data$marker_list$coagulation)
    result_cv_ery <- get_cv(raw_data = data$rawdata,protein = data$marker_list$erythrocyte)
    result_cv_pla <- get_cv(raw_data = data$rawdata,protein = data$marker_list$platelet)
    result_cv_other <- get_cv(raw_data =  data$rawdata[!rownames(data$rawdata)%in%
                                              c(data$marker_list$coagulation,
                                                data$marker_list$erythrocyte,
                                                data$marker_list$platelet),])
    result_cv_coa$type <- "coagulation"
    result_cv_ery$type <- "erythrocyte"
    result_cv_pla$type <- "platelet"
    result_cv_other$type <- "other protein"
    result_cv <- rbind(result_cv_coa,result_cv_ery,result_cv_pla,result_cv_other)
    result_cv$type <- factor(result_cv$type,
                             levels = c("coagulation",
                                        "erythrocyte",
                                        "platelet",
                                        "other protein"))

    ggplot(result_cv,aes(x = type , y = CV, fill = type)) +
      geom_violin() +
      geom_boxplot(fill = "white",width = 0.2) +
      scale_fill_manual(values = c("coagulation" = "#E64B35FF",
                                   "erythrocyte" = "#F39B7FFF",
                                   "platelet" = "#7E6148FF",
                                   "other protein" = "#3C5488FF")) +
      theme_classic() +
      theme(axis.text.x = element_blank())
  })
  output$cv_post_plot <- renderPlot({
    source("./R/modules/get_cv.R")
    req(result_correct())
    data <- result_correct()
    result_cv_coa <- get_cv(raw_data = data$correct_data,protein = data$marker_list$coagulation)
    result_cv_ery <- get_cv(raw_data = data$correct_data,protein = data$marker_list$erythrocyte)
    result_cv_pla <- get_cv(raw_data = data$correct_data,protein = data$marker_list$platelet)
    result_cv_other <- get_cv(raw_data =  data$correct_data[!rownames(data$correct_data)%in%
                                                         c(data$marker_list$coagulation,
                                                           data$marker_list$erythrocyte,
                                                           data$marker_list$platelet),])
    result_cv_coa$type <- "coagulation"
    result_cv_ery$type <- "erythrocyte"
    result_cv_pla$type <- "platelet"
    result_cv_other$type <- "other protein"
    result_cv <- rbind(result_cv_coa,result_cv_ery,result_cv_pla,result_cv_other)
    result_cv$type <- factor(result_cv$type,
                             levels = c("coagulation",
                                        "erythrocyte",
                                        "platelet",
                                        "other protein"))
    
    ggplot(result_cv,aes(x = type , y = CV, fill = type)) +
      geom_violin() +
      geom_boxplot(fill = "white",width = 0.2) +
      scale_fill_manual(values = c("coagulation" = "#E64B35FF",
                                   "erythrocyte" = "#F39B7FFF",
                                   "platelet" = "#7E6148FF",
                                   "other protein" = "#3C5488FF")) +
      theme_classic() +
      theme(axis.text.x = element_blank())
  })
  ### erythrocyte ----
  
  output$contamination_erythrocyte_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_contamination)) {
      return(NULL)  # 如果 plot_contamination 为 NULL，不绘制图表
    }
    plot(result_check()$plot_contamination$erythrocyte)
  })
  ### platelet ----
  output$contamination_platelet_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_contamination)) {
      return(NULL)  # 如果 plot_contamination 为 NULL，不绘制图表
    }
    plot(result_check()$plot_contamination$platelet)
  })
  ### coagulation ----
  output$contamination_coagulation_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_contamination)) {
      return(NULL)  # 如果 plot_contamination 为 NULL，不绘制图表
    }
    plot(result_check()$plot_contamination$coagulation)
  })
  ## 显示标记物可视化 ----
  ### pre ----
  #### erythrocyte marker ----
  output$erythrocyte_marker_pre_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_marker$erythrocyte)) {
      return(NULL)  # 如果 erythrocyte_marker_plot 为 NULL，不绘制图表
    }
    plot(result_check()$plot_marker$erythrocyte)
  })
  
  #### platelet marker ----
  output$platelet_marker_pre_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_marker$platelet)) {
      return(NULL)  # 如果 platelet_marker_plot 为 NULL，不绘制图表
    }
    plot(result_check()$plot_marker$platelet)
  })
  
  #### coagulation marker ----
  output$coagulation_marker_pre_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_marker$coagulation)) {
      return(NULL)  # 如果 coagulation_marker_plot 为 NULL，不绘制图表
    }
    print(result_check()$plot_marker$coagulation)
  })
  ### post(需修改) ----
  #### erythrocyte marker ----
  output$erythrocyte_marker_post_plot <- renderPlot({
    req(result_correct())
    source("./R/plot_protein_by_sample.R")
    plot_protein_by_sample(data = result_correct()$correct_data[rownames(result_correct()$correct_data)%in%result_correct()$marker_list$erythrocyte,])
  })
  
  #### platelet marker ----
  output$platelet_marker_post_plot <- renderPlot({
    req(result_correct())
    source("./R/plot_protein_by_sample.R")
    plot_protein_by_sample(data = result_correct()$correct_data[rownames(result_correct()$correct_data)%in%result_correct()$marker_list$platelet,])
  })
  
  #### coagulation marker ----
  output$coagulation_marker_post_plot <- renderPlot({
    req(result_correct())
    source("./R/plot_protein_by_sample.R")
    plot_protein_by_sample(data = result_correct()$correct_data[rownames(result_correct()$correct_data)%in%result_correct()$marker_list$coagulation,])
  })
  # 显示污染矩阵 ----
  output$data_marker_erythrocyte <- renderDT({
    req(result_check())
    datatable(result_check()$data$erythrocyte, options = list(pageLength = 10))  # 每页显示 10 行
  })
  output$data_marker_coagulation <- renderDT({
    req(result_check())
    datatable(result_check()$data$coagulation, options = list(pageLength = 10))  # 每页显示 10 行
  })
  output$data_marker_platelet <- renderDT({
    req(result_check())
    datatable(result_check()$data$platelet, options = list(pageLength = 10))  # 每页显示 10 行
  })
  

  
  # corrected_plot 显示校正后污染水平可视化 ----
  output$corrected_plot <- renderPlot({
    req(result_correct())
    if (is.null(result_correct()$contamination_level)) {
      return(NULL)  # 如果 contamination_level 为 NULL，不绘制图表
    }
    plot_contamination(result_correct()$contamination_level, "校正后污染水平可视化", "corrected_plot")
  })
  
  # 差异表达分析 ----
  result_de_post <- eventReactive(input$run_de, {
    req(result_correct())
    limma_proteomics_analysis(expr_matrix = log2(result_correct()$correct_data),
                              group_matrix = result_correct()$group,
                              compare = c(input$group1,input$group2),
                              p_type = "raw")
  })
  result_de_pre <- eventReactive(input$run_de, {
    req(result_correct())
    limma_proteomics_analysis(expr_matrix = log2(result_correct()$rawdata),
                              group_matrix = result_correct()$group,
                              compare = c(input$group1,input$group2),
                              p_type = "raw")
  })
  
  # 显示差异表达结果 ----
  output$result_de_pre_table <- renderDT({
    req(result_de_pre())
    datatable(result_de_pre(), options = list(pageLength = 10))  # 每页显示 10 行
  })
  output$result_de_post_table <- renderDT({
    req(result_de_post())
    datatable(result_de_post(), options = list(pageLength = 10))  # 每页显示 10 行
  })
  
  # 绘制火山图 ----
  output$volc_de_pre <- renderPlot({
    req(result_de_pre())
    plot_volc <- create_volcano_plot(result_de_pre(), 
                                     p_type = "raw", 
                                     p_cutoff = 0.05, 
                                     logFC_cutoff = 0,
                                     gene_col = "Protein",
                                     group_names = c(input$group1,input$group2),
                                     colors = c(Up = "#E64B35", Down = "#4DBBD5", Not = "grey80"))
    print(plot_volc)
  })
  output$volc_de_post <- renderPlot({
    req(result_de_post())
    plot_volc <- create_volcano_plot(result_de_post(), 
                                     p_type = "raw", 
                                     p_cutoff = 0.05, 
                                     logFC_cutoff = 0,
                                     gene_col = "Protein",
                                     group_names = c(input$group1,input$group2),
                                     colors = c(Up = "#E64B35", Down = "#4DBBD5", Not = "grey80"))
    print(plot_volc)
  })
  # 下载校正数据 ----
  output$download_data <- downloadHandler(
    filename = function() {
      "corrected_data.csv"
    },
    content = function(file) {
      write.csv(result_correct()$correct_data, file)
    }
  )
  
  # 下载差异表达结果 ----
  ## 矫正前 ----
  output$download_de_pre <- downloadHandler(
    filename = function() {
      "de_results.csv"
    },
    content = function(file) {
      write.csv(result_de_pre(), file)
    }
  )
  ## 矫正后 ----
  output$download_de_post <- downloadHandler(
    filename = function() {
      "de_results.csv"
    },
    content = function(file) {
      write.csv(result_de_post(), file)
    }
  )
}



# 运行 Shiny 应用 ----
shinyApp(ui = ui, server = server)
