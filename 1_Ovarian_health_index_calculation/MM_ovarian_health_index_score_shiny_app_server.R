library(shiny)
library(dplyr)
library(beeswarm)
library(outliers)
library(ggplot2)

YF_medians <- list(
  AMH_median = 133.843,
  FSH_median = 26.02,
  INHBA_median = 405.779,
  follicle_combined_median = 7
)

EF_medians <- list(
  AMH_median = 7.868,
  FSH_median = 23.08,
  INHBA_median = 143.233,
  follicle_combined_median = 0 
)


function(input, output, session) {
  
  ### REFERENCE-BASED CALCULATION
  calculate_scores <- function(value, YF_median, EF_median, increasing) {
    if (increasing) {
      return(case_when(
        value <= YF_median ~ 3,
        value > YF_median & value < EF_median ~ 2,
        TRUE ~ 1
      ))
    } else {
      return(case_when(
        value <= EF_median ~ 1,
        value > EF_median & value < YF_median ~ 2,
        TRUE ~ 3
      ))
    }
  }
  
  process_data <- function(amh, fsh, inhba, follicle) {
    amh_score <- calculate_scores(amh, YF_medians$AMH_median, EF_medians$AMH_median, FALSE)
    fsh_score = calculate_scores(fsh, YF_medians$FSH_median, EF_medians$FSH_median, TRUE)
    inhba_score = calculate_scores(inhba, YF_medians$INHBA_median, EF_medians$INHBA_median, FALSE)
    hormone_avg_score = (amh_score + fsh_score + inhba_score) / 3
    follicle_score = calculate_scores(follicle, YF_medians$follicle_combined_median, EF_medians$follicle_combined_median, FALSE)
    total_score = (hormone_avg_score + follicle_score) / 6 * 100
    return(list(AMH_score=amh_score, FSH_score=fsh_score, INHBA_score=inhba_score, HORMONE_score=hormone_avg_score,
                FOLLICLE_score=follicle_score, TOTAL_score=total_score))
  }

  ### USER-DEFINED CONTROLS
  user_calculate_scores <- function(value, U_YF_median, U_EF_median, increasing) {
    if (increasing) {
      return(case_when(
        value <= YF_median ~ 3,
        value > YF_median & value < EF_median ~ 2,
        TRUE ~ 1
      ))
    } else {
      return(case_when(
        value <= EF_median ~ 1,
        value > EF_median & value < YF_median ~ 2,
        TRUE ~ 3
      ))
    }
  }
  
  user_process_data <- function(amh, fsh, inhba, follicle, U_YF_median, U_EF_median) {
    amh_score <- calculate_scores(amh, U_YF_median$AMH_median, U_EF_median$AMH_median, FALSE)
    fsh_score = calculate_scores(fsh, U_YF_median$FSH_median, U_EF_median$FSH_median, TRUE)
    inhba_score = calculate_scores(inhba, U_YF_median$INHBA_median, U_EF_median$INHBA_median, FALSE)
    hormone_avg_score = (amh_score + fsh_score + inhba_score) / 3
    follicle_score = calculate_scores(follicle, U_YF_median$Follicle_Count_median, U_EF_median$Follicle_Count_median, FALSE)
    total_score = (hormone_avg_score + follicle_score) / 6 * 100
    return(list(AMH_score=amh_score, FSH_score=fsh_score, INHBA_score=inhba_score, HORMONE_score=hormone_avg_score,
                FOLLICLE_score=follicle_score, TOTAL_score=total_score))
  }
  
  results <- reactiveVal(data.frame())
  batch_results <- reactiveVal(data.frame())
  user_batch_results <- reactiveVal(data.frame())
  
  ### SERVER LOGIC FOR CLICK-AND-DRAG CALCULATOR
  observeEvent(input$calculate,
               {
                 req(input$sample_id, input$group_id, input$amh, input$fsh, input$inhba, input$follicle_count)
                 tryCatch({
                   scores <- process_data(input$amh, input$fsh, input$inhba, input$follicle_count)
                   new_row <- data.frame(
                     Sample_ID = input$sample_id,
                     Group_ID = input$group_id,
                     AMH = input$amh,
                     FSH = input$fsh,
                     INHBA = input$inhba,
                     Follicle_Count = input$follicle_count,
                     AMH_Score = scores$AMH_score,
                     FSH_Score = scores$FSH_score,
                     INHBA_Score = scores$INHBA_score,
                     Hormone_Avg_Score = scores$HORMONE_score,
                     Follicle_Score = scores$FOLLICLE_score,
                     Total_Score = scores$TOTAL_score
                   )
                   results(rbind(results(), new_row))
                 }, error = function(e){
                   showNotification(paste('Critical Error! ', e), type = "error")
                 })
               })
  
  observeEvent(input$delete, {
    results(data.frame())
  })
  
  output$score_table <- renderTable(results())
  
  output$score_plot <- renderPlot({
    req(nrow(results()) > 0)
    
    ggplot(results(), aes(y=Total_Score, factor(Group_ID))) +
      geom_boxplot() +
      geom_jitter() +
      labs(x='Group ID', y='Total Score', title='Distribution of Ovarian Health Index')
  })
  
  output$download_results <- downloadHandler(
    filename = function() {
      paste('ovarian_health_index_', Sys.Date(), '.csv')
    },
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
  
  ### SERVER SIDE LOGIC FOR CSV BATCH PROCESSING
  observeEvent(input$process_batch,
               {
                 req(input$file)
                 tryCatch({
                   df <- read.csv(input$file$datapath)
                   result <- df %>%
                     rowwise() %>%
                     mutate(
                       scores = list(process_data(AMH, FSH, INHBA, Follicle_Count)),
                       AMH_score = scores$AMH_score,
                       FSH_score = scores$FSH_score,
                       INHBA_score = scores$INHBA_score,
                       Hormone_Avg_Score = scores$HORMONE_score,
                       Follicle_Score = scores$FOLLICLE_score,
                       Total_Score = scores$TOTAL_score
                     ) %>%
                     select(-scores)
                   batch_results(result)
                 }, error = function(e) {
                   showNotification(paste('Critical Error! ', e), type = "error")
                 }
                 )
               }
  )
  
  output$batch_table <- renderTable(batch_results())
  
  observeEvent(input$batch_delete,
               {
                 batch_results(data.frame())
               })
  
  output$batch_plot <- renderPlot({
    req(nrow(batch_results()) > 0)
    
    ggplot(batch_results(), aes(y=Total_Score, factor(Group_ID))) +
      geom_boxplot() +
      geom_jitter() +
      labs(x='Group ID', y='Total Score', title='Distribution of Ovarian Health Index')
  })
  
  output$download_batch <- downloadHandler(
    filename = function() {
      paste('ovarian_health_index_batch_', Sys.Date(), '.csv')
    },
    content = function(file) {
      write.csv(batch_results(), file, row.names = FALSE)
    }
  )

  ### SERVER SIDE LOGIC FOR USER-DEFINED CONTROL
  observeEvent(input$user_process_batch,
               {
                 req(input$user_file)
                 tryCatch({
                   df <- read.csv(input$user_file$datapath)
                   if(nrow(df %>% filter(Group_ID == "c_YF")) == 0){
                     showNotification("You do not have a valid c_YF group!", type = "error")
                   }
                   if(nrow(df %>% filter(Group_ID == "c_EF")) == 0){
                     showNotification("You do not have a valid c_EF group!", type = "error")
                   }
                   medians <- df %>%
                     group_by(Group_ID) %>%
                     summarise(across(c(AMH, FSH, INHBA, Follicle_Count), median, .names = "{.col}_median"))
                   User_YF_medians <- filter(medians, Group_ID == "c_YF")
                   User_EF_medians <- filter(medians, Group_ID == "c_EF")
                   result <- df %>%
                     rowwise() %>%
                     mutate(
                       scores = list(user_process_data(AMH, FSH, INHBA, Follicle_Count, User_YF_medians, User_EF_medians)),
                       AMH_score = scores$AMH_score,
                       FSH_score = scores$FSH_score,
                       INHBA_score = scores$INHBA_score,
                       Hormone_Avg_Score = scores$HORMONE_score,
                       Follicle_Score = scores$FOLLICLE_score,
                       Total_Score = scores$TOTAL_score
                     ) %>%
                     select(-scores)
                   user_batch_results(result)
                 }, error = function(e){
                   showNotification(paste('Critical Error! ', e), type = "error")
                 }
                 )
               }
  )
  
  observeEvent(input$user_batch_delete,
               {
                 user_batch_results(data.frame())
               })
  
  output$user_batch_table <- renderTable(user_batch_results())
  
  output$user_batch_plot <- renderPlot({
    req(nrow(user_batch_results()) > 0)
    
    ggplot(user_batch_results(), aes(y=Total_Score, factor(Group_ID))) +
      geom_boxplot() +
      geom_jitter() +
      labs(x='Group ID', y='Total Score', title='Distribution of Ovarian Health Index')
  })
  
  output$user_download_batch <- downloadHandler(
    filename = function() {
      paste('ovarian_health_index_user_batch_', Sys.Date(), '.csv')
    },
    content = function(file) {
      write.csv(user_batch_results(), file, row.names = FALSE)
    }
  )
}


