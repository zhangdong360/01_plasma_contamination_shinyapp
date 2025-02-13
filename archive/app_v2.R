library(shiny)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(readxl)
library(DT)  # 引入 DT 包
source("./data_check_v02.R")

# 定义 UI
ui <- fluidPage(
  titlePanel("血浆蛋白污染校正与差异表达分析"),
  sidebarLayout(
    sidebarPanel(
      fileInput("data_file", "上传数据文件 (CSV)", accept = ".csv"),
      fileInput("group_file", "上传分组信息文件 (CSV)", accept = ".csv"),
      selectInput("type", "校正类型", choices = c("all", "erythrocyte", "platelet", "coagulation")),
      selectInput("group1", "分组1", choices = NULL),
      selectInput("group2", "分组2", choices = NULL),
      actionButton("run_check", "运行数据检查"),
      actionButton("run_correct", "运行数据校正"),
      actionButton("run_de", "运行差异表达分析"),
      downloadButton("download_data", "下载校正数据"),
      downloadButton("download_de", "下载差异表达结果")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("原始数据", 
                 h3("数据文件"),
                 DTOutput("data_table"),  # 使用 DTOutput 替换 tableOutput
                 h3("分组信息文件"),
                 DTOutput("group_table")),  # 使用 DTOutput 替换 tableOutput
        tabPanel("数据检查结果", 
                 h3("缺失基因"),
                 verbatimTextOutput("missing_genes"),
                 # 二级选项卡：污染情况
                 h3("污染marker表达情况"),
                 tabsetPanel(
                   tabPanel("红系污染",
                            h4("matrix"),
                            DTOutput("data_marker_erythrocyte"),
                            h4("marker plot"),
                            plotOutput("contamination_erythrocyte_plot")),
                   tabPanel("凝血污染", 
                            h4(" matrix"),
                            DTOutput("data_marker_coagulation"),
                            h4("marker plot"),
                            plotOutput("contamination_coagulation_plot")),
                   tabPanel("血小板污染", 
                            h4("matrix"),
                            DTOutput("data_marker_platelet"),
                            h4("marker plot"),
                            plotOutput("contamination_platelet_plot"))
                 ),
                 h3("样本污染情况"),
                 tabsetPanel(
                   tabPanel("红系污染",
                            plotOutput("erythrocyte_marker_plot")),
                   tabPanel("凝血污染", 
                            plotOutput("coagulation_marker_plot")),
                   tabPanel("血小板污染", 
                            plotOutput("platelet_marker_plot"))
                 )
        ),
        tabPanel("校正结果", 
                 h3("校正后污染水平可视化"),
                 plotOutput("corrected_plot")),
        tabPanel("差异表达分析", 
                 h3("差异表达结果"),
                 DTOutput("de_table"))  # 使用 DTOutput 替换 tableOutput
      )
    )
  )
)

# 定义 Server 逻辑
server <- function(input, output, session) {
  # 读取数据
  data <- reactive({
    req(input$data_file)
    df <- read.csv(input$data_file$datapath)
    rownames(df) <- df[, 1]  # 将第一列设置为行名
    df <- df[, -1, drop = FALSE]  # 删除第一列
    df
  })
  
  data_group <- reactive({
    req(input$group_file)
    df <- read.csv(input$group_file$datapath)
    rownames(df) <- df[, 1]  # 将第一列设置为行名
    df
  })
  
  # 展示原始数据（使用 DT 包）
  output$data_table <- renderDT({
    req(data())
    datatable(data(), options = list(pageLength = 10))  # 每页显示 10 行
  })
  
  output$group_table <- renderDT({
    req(data_group())
    datatable(data_group(), options = list(pageLength = 10))  # 每页显示 10 行
  })
  
  # 更新分组选择
  observe({
    req(data_group())
    updateSelectInput(session, "group1", choices = unique(data_group()$group))
    updateSelectInput(session, "group2", choices = unique(data_group()$group))
  })
  
  # 数据检查
  result_check <- eventReactive(input$run_check, {
    data_check(data(), data_group())
  })
  
  # 显示缺失基因
  output$missing_genes <- renderPrint({
    req(result_check())
    cat("红系marker缺失基因:\n")
    print(result_check()$missing_genes$missing_erythrocyte_genes)
    cat("凝血marker缺失基因:\n")
    print(result_check()$missing_genes$missing_coagulation_genes)
    cat("血小板marker缺失基因:\n")
    print(result_check()$missing_genes$missing_platelet_genes)
  })
  
  # 显示污染水平可视化
  output$contamination_erythrocyte_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_contamination)) {
      return(NULL)  # 如果 plot_contamination 为 NULL，不绘制图表
    }
    # 绘制污染水平可视化
    plot(result_check()$plot_contamination$erythrocyte)
  })
  output$contamination_platelet_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_contamination)) {
      return(NULL)  # 如果 plot_contamination 为 NULL，不绘制图表
    }
    # 绘制污染水平可视化
    plot(result_check()$plot_contamination$platelet)
  })
  output$contamination_coagulation_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_contamination)) {
      return(NULL)  # 如果 plot_contamination 为 NULL，不绘制图表
    }
    # 绘制污染水平可视化
    plot(result_check()$plot_contamination$coagulation)
  })
  
  # 显示红细胞标记物可视化
  output$erythrocyte_marker_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_marker$erythrocyte)) {
      return(NULL)  # 如果 erythrocyte_marker_plot 为 NULL，不绘制图表
    }
    plot(result_check()$plot_marker$erythrocyte)
  })
  
  # 显示血小板标记物可视化
  output$platelet_marker_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_marker$platelet)) {
      return(NULL)  # 如果 platelet_marker_plot 为 NULL，不绘制图表
    }
    plot(result_check()$plot_marker$platelet)
  })
  
  # 显示凝血标记物可视化
  output$coagulation_marker_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_marker$coagulation)) {
      return(NULL)  # 如果 coagulation_marker_plot 为 NULL，不绘制图表
    }
    print(result_check()$plot_marker$coagulation)
  })
  
  # 显示污染矩阵
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
  
  # 数据校正
  result_correct <- eventReactive(input$run_correct, {
    data_correct(result_check(), type = input$type)
  })
  
  # 显示校正后污染水平可视化
  output$corrected_plot <- renderPlot({
    req(result_correct())
    if (is.null(result_correct()$contamination_level)) {
      return(NULL)  # 如果 contamination_level 为 NULL，不绘制图表
    }
    plot_contamination(result_correct()$contamination_level, "校正后污染水平可视化", "corrected_plot")
  })
  
  # 差异表达分析
  result_de <- eventReactive(input$run_de, {
    DE_analysis(result_correct(), type = "post", group_1 = input$group1, group_2 = input$group2)
  })
  
  # 显示差异表达结果
  output$de_table <- renderDT({
    req(result_de())
    datatable(result_de(), options = list(pageLength = 10))  # 每页显示 10 行
  })
  
  # 下载校正数据
  output$download_data <- downloadHandler(
    filename = function() {
      "corrected_data.csv"
    },
    content = function(file) {
      write.csv(result_correct()$correct_data, file)
    }
  )
  
  # 下载差异表达结果
  output$download_de <- downloadHandler(
    filename = function() {
      "de_results.csv"
    },
    content = function(file) {
      write.csv(result_de(), file)
    }
  )
}

# 运行 Shiny 应用
shinyApp(ui = ui, server = server)
