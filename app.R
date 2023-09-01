## Application intends to generate a molecular profile for the binary classification of samples from a processed omics dataset
# Data should be in the format described in the first tab

# Libraries -----
library(shiny)
library(DT)
library(vroom)
library(caret)
library(ggplot2)
library(parallel)
library(doParallel)
library(dplyr)
library(plyr)
library(magrittr)
library(data.table)
library(MLeval)
library(limma)
library(purrr)
library(kableExtra)
library(ggpubr)
library (gtools)
library(tidyverse) 
library(ggrepel)
library(tidyr)
library(e1071)
library(randomForest)
library(ggplot2)
library(tidyverse) 
library (readr)
library(PMCMRplus)
library (PMCMR)
library(ggpubr)
library (mosaic)
library (dplyr)
library (data.table)
library(reshape2)
library (gtools)
library(plyr)
library(limma)
library(ggrepel)
library(amap)
library(rstatix)
library(broom)
library(ggprism)
library(HDDesign)
library(caret)
library(rsample)
library(sandwich)
library(rpart)
library(rpart.plot)
library(randomForest)
library(RColorBrewer)
library(plotly)
library(purrr)
library(e1071)
library(ggraph)
library(igraph)
library(pscl)
library(parallel)
library(doParallel)
library(ROCR) 
library(corrr)
library(ggcorrplot)
library(ape)
library(forcats)
library(kernlab)
library(xgboost)
library(mlbench)
library(naivebayes)
library(png)
library(DALEX)
library(DALEXtra)



#Helpful functions ------
`%notin%` <- Negate(`%in%`)
`%notlike%` <- Negate(`%like%`)
back_2_front <- function(df){
  df <- df[,c(ncol(df),1:(ncol(df)-1))]
  return(df)
}

flextable_only <- function(table){
  table_flex <- flextable::flextable(table)%>%
    flextable::fontsize(size = 8)%>%
    flextable::autofit()%>%
    flextable:: height( height = .5)%>%
    flextable::theme_vanilla()
}

# deploy app to github using: rsconnect::deployApp('ML_taser/')
options(repos = BiocManager::repositories()) 

library(rsconnect)
rsconnect::deployApp()

thematic::thematic_shiny(font = "auto")

# Define UI -----
ui <- fluidPage( 
  navbarPage(
    theme = shinythemes::shinytheme("united"),
    HTML('<a style="text-decoration:none; 
         cursor:default;
         color:#FFFFFF;
         font-size:30px; 
         font-family:arial; 
         text-align: right;" 
         class="active" 
         href="#">School of Infection and Immunity, University of Glasgow</a>'),
    tags$script(HTML("var header = $('.navbar > .container-fluid'); header.append('<a href=\"URL\"><img src=\"UofG_Coat_of_Arms.png\"  style=\"float:right;width:100px;height:150px;padding-top:10px;\"> </a></div>');console.log(header)")),
  
    id="nav"),

  titlePanel(div(column(width = 10, 
                       h2("Generation of a Molecular Profile of Binary Classes from an Omic Dataset")))),
  
  mainPanel(div(column(width = 12, 
                       tabsetPanel(
                   tabPanel("Import File",
                            mainPanel(
                              column(width=5),
                              # Importing file
                              tags$hr(),
                              h2('Import Omic Data File'),
                              h5("See 'Data and Formatting' tab for example dataset and formatting"),
                              tags$hr(),
                              fileInput("upload", NULL, accept = c(".csv", ".tsv")),
                              h2("Study Groups in Dataset"),
                              h4("Define binary groups for analysis and molecular profile generation"),
                              h5("e.g. study condition 1: healthy, study condition 2: disease"),
                              textInput("condition_1", 'Study Condition 1'),
                              textOutput('cond_1'),
                              tags$hr(),
                              textInput("condition_2", 'Study Condition 2'),
                              textOutput('cond_2'),
                              tags$hr(),
                              h2("Overview of Imported Dataset"),
                              numericInput("n", "Select n samples to show", value = 5, min = 1, step = 1),
                              tags$hr(),
                              tableOutput("head")
                              )),
                   
                   tabPanel("Exploratory Data Analysis",
                            tags$hr(),
                            h2("Condition Numbers in Dataframe"),
                            plotOutput("class_balance_actual"),
                            tags$hr(),
                            h2("PCA plot"),
                            plotOutput("PCA"),
                            tags$hr(),
                            h2("Volcano Plot: Differential Abundance of Features Across Study Conditions"),
                            plotOutput("limma"),
                            tags$hr(),
                            h2("Differential Abundance Table: Limma Output"),
                            tableOutput("limma_tab")
                   ),
                   tabPanel("Model Training",
                            verticalLayout(
                              sliderInput("train_test_split", "Select train:test subset split (e.g. select 0.7 for 70:30 training: testing split)", min= 0.5, max= 1.0, value = 0.7, step= 0.05),
                              numericInput('mod_gen_repeats', "Enter number of repeats K-fold cross validations used in resampling: ", value= 1, min=1, step=1),
                              numericInput('repeatedcv_k', "Number of folds in K-fold cross-validation (e.g. 50-100)", value= 10, min=1, step=1),
                              numericInput('fs', 'Number of features used', value= 10, min=1, step=1),
                              tags$hr(),
                              h3('Please be patient, model training can take time...'),
                              tags$hr(),
                              h2("Feature Selection"),
                              plotOutput("gg_fs"),
                              tags$hr(),
                              h2("Algorithm Selection"),
                              fluidRow(
                                splitLayout(cellWidths = c("50%", "50%"), 
                                            plotOutput("model_selection_roc"),
                                            plotOutput("model_selection_pr")))),
                            tableOutput('algorithm_selection_table'),
                            tags$hr(),
                            h2("Performance Metrics"),
                            selectInput('model_from_train', 'Select ML algorithm', c('svmRadial', 'rf', 'xgbTree', 'knn', 'naive_bayes', 'glm', 'glmboost')),
                            plotOutput("model_train_roc")),
                   
                   tabPanel("Model Evaluation",
                          mainPanel(
                              h1("Evaluation of Final Model"),
                              tags$hr(),
                              #tableOutput("head_test"),
                              verticalLayout(
                                h2("ROC Curve"),
                                tags$hr(),
                                plotOutput("model_test_roc"),
                                h2("Calibration Curve"),
                                plotOutput("model_test_cc"),
                                tags$hr(),
                                h2('Final model metrics'),
                                tableOutput("Test_table")))),
                   
                   tabPanel("Interpretting Model Features",
                            h1("Influence of the Selected Features in Model"),
                            tags$hr(),
                            h2('Accumulated Local Effects Plots and Comparison with Partial Dependence Plots'),
                            tags$hr(),
                            plotOutput("ft_int_ale_plot"),
                            tags$hr(),
                            # h2("Shapley Additive Explanations Plot of Selected Features in Model"),
                            # tags$hr(),
                            # plotOutput("shapp"),
                            # tags$hr(),
                            h2("Actual Levels of the Features in Model Across Sample Groups"),
                            plotOutput("actual_plots"),
                            tags$hr(),
                            tags$hr(),
                            tags$hr(),
                            plotOutput('combined_plots')
                            ),
                   
                   tabPanel("Data and Formatting",
                            h2("Example Dataset"),
                            tags$hr(),
                            h3('First 10 Features'),
                            
                            mainPanel(
                              tableOutput("example_table"),
                              tags$hr(),
                              h2("Class Balance: Study Conditions"),
                              plotOutput("class_balance"),
                              tags$hr()
                            ))
                 )))))

# Define server ----- 
server <- function(input, output, session) ({
  thematic::thematic_shiny()
  load_file <- function(name, path) {
    ext <- tools::file_ext(name)
    switch(ext,
           csv = example_df,
           #tsv = vroom::vroom(path, delim = "\t"),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
  }
  
  example_df <- reactive({
    example_df <- data.frame(ID = paste0('Feature_', 1:200),
                             control_1= sample(0:1000, size=200, replace=T),
                             control_2= sample(0:1200, size=200, replace=T),
                             control_3= sample(0:1800, size=200, replace=T),
                             control_4= sample(0:1000, size=200, replace=T),
                             control_5= sample(0:1200, size=200, replace=T),
                             control_6= sample(0:1800, size=200, replace=T),
                             control_7= sample(0:1000, size=200, replace=T),
                             control_8= sample(0:1200, size=200, replace=T),
                             control_9= sample(0:1800, size=200, replace=T),
                             disease1_1= sample(700:3300, size=200, replace=T), 
                             disease1_2= sample(1000:2000, size=200, replace=T),
                             disease1_3= sample(900:2200, size=200, replace=T),
                             disease1_4 = sample(0:2000, size=200, replace=T),
                             disease2_1 = sample(0:1800, size=200, replace=T),
                             disease2_2 = sample(0:1000, size=200, replace=T),
                             disease2_3 = sample(0:2000, size=200, replace=T),
                             disease2_4 = sample(0:2500, size=200, replace=T),
                             disease2_5 = sample(0:1500, size=200, replace=T))
    rownames(example_df) <- example_df$ID
    example_df <- example_df[,-1]
    example_df <- as.data.frame(t(example_df))
    example_df$Condition <- rownames(example_df)
    example_df <- example_df %>% 
      relocate(Condition)
    example_df$Condition <- gsub("(.*)_.*","\\1",example_df$Condition)
    example_df$Condition <- as.factor(example_df$Condition)
    example_df
    
  })
  
  length_example_df <- reactive({
    length(example_df()$Condition)
  })
  
  output$example_table <- renderTable({
    example_df()[1:length_example_df(), 1:11]
  })
  
  data <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ","),
           tsv = vroom::vroom(input$upload$datapath, delim = "\t"),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  output$condition_names <- renderText({
    condition_names <- as.factor(data()$Condition)
    conditions <- levels(condition_names)
    conditions
  })
  
  
  output$head <- renderTable({
    req(input$upload)
    
    head(data()[1:10], input$n)
  })
  
  condition_one <- reactive({
    cond_1_df <- data$Condition %like% input$condition_1
    one <- levels(as.factor(cond_1_df$Condition))
    one
  })

  condition_two <- reactive({
    cond_2_df <- data$Condition %like% input$condition_2
    two <- levels(as.factor(cond_2_df$Condition))
    two
  })
  
  data_use <- reactive({
    data <- data()
    # data <- data %>% 
    #   filter(Condition %like% condition_one() | Condition %like% condition_two())
    data
  })
  
  data_df <- reactive({
    data_df <- data()
    data_df$Condition <- as.factor(data_df$Condition)
    data_df <- data_df[,-2]
    data_df
  })
  
  condition_1_check <- reactive({
    conditions <- levels(data_df()$Condition)
    
    cond_check <- ifelse(conditions[1] == input$condition_1, 'Condition 1 Matches Dataset', "Make sure first condition name matches that described in dataset")
    cond_check
  })
  
  output$cond_1 <- renderText({
    condition_1_check()
  })
  
  condition_2_check <- reactive({
    conditions <- levels(data_df()$Condition)
    
    cond_check <- ifelse(conditions[2] == input$condition_2, 'Condition 2 Matches Dataset', "Make sure second condition name matches that described in dataset")
    cond_check
  })
  
  output$cond_2 <- renderText({
    condition_2_check()
  })
  
  data_df_eg <- reactive({
    data_df <- data_use()
    data_df$Condition <- as.factor(data_df$Condition)
    data_df
  })
  
  output$head_data <- renderTable({
    data_df_eg()[1:10, 1:11]
  })
  
  data_train <- reactive({
    set.seed(42)
    index <- createDataPartition(data_df()$Condition, p = input$train_test_split, list = FALSE)
    train_data <- data_df()[index, ]
  })
  
  data_test <- reactive({
    set.seed(42)
    index <- createDataPartition(data_df()$Condition, p = input$train_test_split, list = FALSE)
    test_data <- data_df()[-index,]
  })
  
  output$class_balance_actual <- renderPlot({
    ggplot(data_df())+
      geom_bar(aes(x=Condition), fill= c('#E69F00', '#56B4E9'))+
      theme_minimal() +
      theme(axis.text= element_text(size=14),
            axis.title=element_text(size=16))
  })
  
  PCA_data <- reactive({
    comb_ft_PCA <- data_df()[,-c(1)]
    scaled_intensities <- scale((comb_ft_PCA))
    scaled_intensities[do.call(cbind, lapply(scaled_intensities, is.nan))] <- 0
    scaled_intensities<- as.data.frame(scaled_intensities)
    pca_data <- prcomp(scaled_intensities)
    pca_data
  })
  
  output$PCA <- renderPlot({
    pca_coord <- data.frame(PCA_data()$x)
    var_explained <- PCA_data()$sdev^2/sum(PCA_data()$sdev^2)
    var_explained[1:5]
    pca_coord$Condition <-as.factor(data_df()$Condition)
    ggplot(pca_coord) + 
      geom_point(size=5, alpha=0.7, 
                 aes(x=PC1,y=PC2, colour= Condition, fill= Condition))+
      labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
           y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))+
      geom_hline(yintercept = 0,
                 colour='navy',
                 linetype='dashed')+
      geom_vline(xintercept = 0,
                 colour='navy',
                 linetype='dashed')+
      theme_minimal() +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=16),
            legend.text = element_text(size=16),
            legend.title = element_text(size=16))
  })

  output$class_balance <- renderPlot({
    ggplot(example_df())+
      geom_bar(aes(x=Condition), fill= c('#E69F00', '#56B4E9', 'chartreuse1'))+
      theme_minimal() +
      theme(axis.text= element_text(size=14),
            axis.title=element_text(size=16))
  })
  
  
  # Differential abundance ----
  volcano_df <- reactive({
    volc <- data_df()[,-c(1)]
    volc$Condition <- as.factor(data_df()$Condition)
    volc <- back_2_front(volc)
    rownames(volc) <- paste0(volc$Condition, '_', rownames(volc))
    volc
  })
  
  limma_fun <- reactive({
    df1 <- volcano_df()
    df1$Condition <- ifelse(df1$Condition == input$condition_1, 'no', 'yes')
    conditions <- levels(as.factor(df1$Condition))
    
    rownames(df1) <- paste0(df1$Condition, '_', c(1:length(df1$Condition)))
  
    df2 <- df1[,-1]
    df2 <- as.data.frame(t(df2))
    colnames(df2) <-  df1$Condition
    df2[1:10,1:10]
    
    Group <- factor(names(df2), levels = c('no', 'yes'))
    design <- model.matrix (~Group)
    colnames(design) <- c('no', 'novsyes')
    eset <- df2
    fit <- lmFit(eset, design)
    fit <- eBayes(fit)
    toptable <- topTable(fit, coef = 'novsyes', adjust = 'BH', number = 1800)
    toptable <- as.data.frame(toptable)
    toptable$Feature <- rownames(toptable)
    toptable <- toptable[,c(ncol(toptable),1:(ncol(toptable)-1))]
    toptable$Sig <- 0
    toptable$Sig <- ifelse(toptable$adj.P.Val <0.05, '< 0.05', '> 0.05')
    toptable$Sig_Names <-0
    toptable$Sig_Names <- ifelse(toptable$Sig =='< 0.05' ,toptable$Feature, '')
    toptable
    
  })
  
  output$limma_tab <- renderTable({
    head(limma_fun()[c(1:7)], n=10)
  })
  
  output$limma <- renderPlot({
    ggplot(limma_fun(), aes(x=logFC, y=-log10(P.Value), 
                            colour=Sig, 
                            group=Sig)) +
      geom_point (size=3, alpha=0.7) +
      theme_minimal() +
      labs (x='LogFC',
            y='-Log p-value',
            colour='Adjusted \np-value')+
      # geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Names),
      #                 show.legend = F,
      #                 box.padding =1,
      #                 size=5,
      #                 max.overlaps = Inf,
      #                 position = position_jitter(seed = 1),
      #                 arrow = arrow(length = unit(0.0015, "npc"))) +  
      theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            legend.text= element_text(size=16),
            legend.title= element_text(size=14)) +
      scale_color_brewer(palette = "Set1",direction=-1)       
    })
  
  #Machine learning ----
  ft_sel <- reactive({
    options(warn=-1)
    subsets <- c(1:10)
    set.seed(42)   
    ctrl <- rfeControl(functions = rfFuncs,
                       method = "repeatedcv",
                       number = 10,
                       repeats=1,
                       verbose = FALSE)
    
    profile <- rfe(x=data_train()[,-1], y=data_train()$Condition,
                   sizes = subsets,
                   rfeControl = ctrl)
    profile_df <- profile$variables
    profile_df[1:input$fs,]
      })
  
  output$gg_fs <- renderPlot({
    ab <- ft_sel()
    ggplot(ab, aes(x=Overall, y=reorder(var, Overall)))+
      geom_col(aes(fill=Overall))+
      theme_minimal()+
      scale_fill_continuous(low='light green', high='navy')+
      theme(legend.position = 'none',
            axis.title.y = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10))+
      labs(x='Relative Importance', 
           y='Feature')
  })
  
  model_df <- reactive({
    df <- data_df() %>% 
      dplyr::select(Condition, ft_sel()$var)
      df
  })
  
  model_sel <- reactive({
    control= trainControl(method="repeatedcv", 
                          number=input$repeatedcv_k, 
                          repeats=input$mod_gen_repeats,
                          summaryFunction=twoClassSummary,
                          savePredictions = TRUE, 
                          classProbs = TRUE, 
                          verboseIter = TRUE,
                          search = 'random')
    
    tunelength= 15
    data_df = model_df()
    # train the linear model
    set.seed(42)
    modelLm <- caret::train(Condition~., data=data_df, method='glm', trControl=control, tuneLength=tunelength)
    # train the KNN model
    set.seed(42)
    modelKNN <- caret::train(Condition~., data=data_df, method='kknn', trControl=control, tuneLength=tunelength)
    # train Boosted GLM model
    set.seed(42)
    modelGLM <- caret::train(Condition ~., data=data_df, method='glmboost',trControl=control,tuneLength=tunelength)
    # train the SVM model
    set.seed(42)
    modelSvm <- caret::train(Condition~., data=data_df, method="svmRadial", trControl=control,tuneLength = tunelength)
    # train the RF model
    set.seed(42)
    modelRF <- caret::train(Condition ~. ,data=data_df, method="ranger", trControl=control, tuneLength = tunelength)
    # train the XGB model
    set.seed(42)
    modelXGB <- caret::train(Condition~., data=data_df, method='xgbTree',trControl=control,tuneLength=tunelength)
    
    comp_roc <- MLeval::evalm(list(modelLm, modelKNN, modelGLM, modelSvm,  modelRF, modelXGB),
                      gnames=c('LM', 'KNN', 'GLM', 'SVM', 'RF', "XGB"))
    comp_roc
  }) 
  
  alg_sel_table <- reactive({
    ml_eval_output <- as.data.frame(model_sel()$stdres)
    ml_eval_output$Measure <- rownames(ml_eval_output)
    ml_eval_output <- back_2_front(ml_eval_output)
    ml_eval_output <- flextable_only(ml_eval_output)
  })
  
  output$model_selection_roc <- renderPlot({ # plot ROC for all models from training data
    model_sel()[[1]]
  })
  
  output$model_selection_pr <- renderPlot({ # plot ROC for all models from training data
    model_sel()[[2]]
  })
  
  alg_sel_tb <- reactive({
    ml_eval_output <- as.data.frame(model_sel()$stdres)
    ml_eval_output$Measure <- rownames(ml_eval_output)
    ml_eval_output <- back_2_front(ml_eval_output)
    ml_eval_output[c(1:3,8:14),]
  })
  
  output$algorithm_selection_table <- renderTable({
    alg_sel_tb()
  })
  
  ## Model generation
  model_gen <- reactive({
    set.seed(42)
    modell <- train(Condition~., data=model_df(),
                   method=input$model_from_train, 
                    metric='ROC',
                    trControl= trainControl(method="repeatedcv", 
                                            number=10, 
                                            repeats=10,
                                            summaryFunction=twoClassSummary,
                                            savePredictions = TRUE, 
                                            classProbs = TRUE, 
                                            verboseIter = TRUE))
    modell
  })
  
  output$model_train_roc <- renderPlot({
    evalm(model_gen())$roc
  })
  
  test_df_values <- reactive({
      predictions <- predict(model_gen(), data_test())
      con_matr <- confusionMatrix(predictions, data_test()$Condition)
      con_stats <- con_matr$overall
      
      pr <- prediction(as.numeric(predictions), as.numeric(data_test()$Condition))
      prf <- performance(pr, measure = "tpr", x.measure = "fpr")
      auc <- performance(pr, measure = "auc")
      auc_val <- auc@y.values[[1]]
      result.predicted.prob <- predict(model_gen(), data_test(), type="prob") # Prediction
      result.roc <- roc(data_test()$Condition, result.predicted.prob$Control) # Draw ROC curve.
      result.roc
  })
  
  test_df_actual <-  reactive({
    roc_met_df <- data.frame(test_df_values()$sensitivities,test_df_values()$specificities)
    names(roc_met_df) <- c("Sensitivity", 'Specificity')
    roc_met_df$One_Minus_Spec <- 1-roc_met_df$Specificity
    roc_met_df$Model <- 'Omic Model' 
    roc_met_df
  })
  
  output$model_test_roc <- renderPlot({ 
    pred <- predict(model_gen(), newdata=data_test(), type="prob")
    
    evalm(data.frame(pred, data_test()$Condition))$roc
    })
  
 
  output$model_test_cc <- renderPlot({ 
    pred <- predict(model_gen(), newdata=data_test(), type="prob")
    
    evalm(data.frame(pred, data_test()$Condition))$cc
  })   
  
  test_table <- reactive({
    pred <- predict(model_gen(), newdata=data_test(), type="prob")
    test_evalm <- evalm(data.frame(pred, data_test()$Condition))
    test_table <- as.data.frame(test_evalm$optres)
    test_table$Metric <- rownames(test_table)
    test_table <- back_2_front(test_table)
    names(test_table) <- c('Metric', 'Score', 'CI')
    test_table
  })
  
  output$Test_table <- renderTable({
    test_table()
    
  })

# Feature Interpretation
  ft_int_behind_scenes <- reactive({
    library(DALEXtra)
    train_data_3 <-  model_df()
    xplainer_rf <- DALEX::explain(model = model_gen(),  
                                  data = train_data_3[, -1],
                                  y = train_data_3$Condition, 
                                  type='classification')
    xplainer_rf
  })
  
  output$ft_int_ale_plot <- renderPlot({
    train_data_3 <- model_df()
    explainer <- ft_int_behind_scenes()
    pd_rf <- DALEX::model_profile(explainer = explainer,
                           type = "partial",
                           variables = names(train_data_3)[2:ncol(train_data_3)])
    ale_rf <- DALEX::model_profile(explainer = explainer, 
                            type='accumulated',
                            variables = names(train_data_3)[2:ncol(train_data_3)])
    
     # the B=25 is the number of random orderings of the explanatory variables
    
    #shapleys_rf <- predict_parts(explainer = xplainer_rf, type='shap',new_observation =  train_data_3, B=25) # the B=25 is the number of random orderings of the explanatory variables
    
    pd_rf$agr_profiles$`_label_` = "partial dependence"
    ale_rf$agr_profiles$`_label_` = "accumulated local"
    
    partials <- plot(pd_rf,
                     ale_rf)
    
    partials
  })
  
  output$shapp <- renderPlot({
    train_data_3 <- model_df()
    explainer <- ft_int_behind_scenes()
    shapleys_rf <- predict_parts(explainer = explainer, 
                                 type='shap',
                                 new_observation =  train_data_3, 
                                 B=25)
    plot(shapleys_rf)
  })

  output$actual_plots <- renderPlot({
    model_df() %>% 
      pivot_longer(cols=2:ncol(model_df()),
                   names_to='Feature',
                   values_to='Level') %>% 
      ggplot(aes(x=Condition, y=Level, fill=Condition)) +
      geom_boxplot()+
      geom_jitter()+
      facet_wrap(~Feature, scale="free_y")+
      stat_compare_means(vjust=2)+
      theme_minimal()+
      theme(axis.text=element_text(size=14),
            strip.text = element_text(size=20),
            axis.title=element_text(size=16),
            legend.position='none')
  })
  
output$combined_plots <- renderPlot({
  # classes
  class_plot <- ggplot(data_df())+
    geom_bar(aes(x=Condition), fill= c('#E69F00', '#56B4E9'))+
    theme_minimal() +
    theme(axis.text= element_text(size=14),
          axis.title=element_text(size=16))
  
  # volcano plot
  volc_plot <- ggplot(limma_fun(), aes(x=logFC, y=-log10(P.Value), 
                          colour=Sig, 
                          group=Sig)) +
    geom_point (size=3, alpha=0.7) +
    theme_minimal() +
    labs (x='LogFC',
          y='-Log p-value',
          colour='Adjusted \np-value')+
    # geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Names),
    #                 show.legend = F,
    #                 box.padding =1,
    #                 size=5,
    #                 max.overlaps = Inf,
    #                 position = position_jitter(seed = 1),
    #                 arrow = arrow(length = unit(0.0015, "npc"))) +  
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.text= element_text(size=16),
          legend.title= element_text(size=14)) +
    scale_color_brewer(palette = "Set1",direction=-1)       
  
  # pca
  pca_plot <- pca_coord <- data.frame(PCA_data()$x)
  var_explained <- PCA_data()$sdev^2/sum(PCA_data()$sdev^2)
  var_explained[1:5]
  pca_coord$Condition <-as.factor(data_df()$Condition)
  pca_plot <- ggplot(pca_coord) + 
    geom_point(size=5, alpha=0.7, 
               aes(x=PC1,y=PC2, colour= Condition, fill= Condition))+
    labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
         y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))+
    geom_hline(yintercept = 0,
               colour='navy',
               linetype='dashed')+
    geom_vline(xintercept = 0,
               colour='navy',
               linetype='dashed')+
    theme_minimal() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16))
  
  # feature selection
  ab <- ft_sel()
  ft_sel <- ggplot(ab, aes(x=Overall, y=reorder(var, Overall)))+
    geom_col(aes(fill=Overall))+
    theme_minimal()+
    scale_fill_continuous(low='light green', high='navy')+
    theme(legend.position = 'none',
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))+
    labs(x='Relative Importance', 
         y='Feature')
  
  # model performance
  pred <- predict(model_gen(), newdata=data_test(), type="prob")
  
  evalm(data.frame(pred, data_test()$Condition))$roc
  
  # feature interpretation
  train_data_3 <- model_df()
  explainer <- ft_int_behind_scenes()
  pd_rf <- DALEX::model_profile(explainer = explainer,
                                type = "partial",
                                variables = names(train_data_3)[2:ncol(train_data_3)])
  ale_rf <- DALEX::model_profile(explainer = explainer, 
                                 type='accumulated',
                                 variables = names(train_data_3)[2:ncol(train_data_3)])
  
  # the B=25 is the number of random orderings of the explanatory variables
  
  #shapleys_rf <- predict_parts(explainer = xplainer_rf, type='shap',new_observation =  train_data_3, B=25) # the B=25 is the number of random orderings of the explanatory variables
  
  pd_rf$agr_profiles$`_label_` = "partial dependence"
  ale_rf$agr_profiles$`_label_` = "accumulated local"
  
  partials <- plot(pd_rf,
                   ale_rf)
  
  partials
  
  # combine all
  library(patchwork)
  class_plot/ 
    pca_plot/ 
    volc_plot/
    ft_sel
})
  
})

shinyApp(ui = ui, server = server)

