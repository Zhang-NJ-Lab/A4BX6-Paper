rm(list = ls()) 
setwd("~/NUIST/John/胡文广/A4BX6/") 
source('JAP/code/dependencies.R')
cl <- mlr_learners$keys()[str_detect(mlr_learners$keys(), 'classif')]

#==== Utilities for this specific task ====
clean_up <- function(df, abbs){
  # df = st1
  df <- df %>%
    select(-contains('Symbol'), -contains('Phase'), -contains('Radioactive'), 
           -contains('Natural'), -contains('Metal'), -contains('Nonmetal'), 
           -contains('Metalloid'), -contains('Type'),-contains('Discoverer'),
           -contains('Group'), -contains('Year'), -contains('NumberofValence'))
  df <- df[!is.na(df$homo), ]
  colnames(df) <- abbs$Abbreviation
  df <- as.data.frame(df)
  st <- df[, -1]
  for(s in 1:ncol(st)){
    st[, s] <- as.numeric(st[, s])
  }
  return(st)
}
get_tune_expression <- function(learner){
  paras_to_tune <- learner$param_set$ids()[!(learner$param_set$class %in% c("Paramlgl", "ParamUty", "ParamFct"))]
  ps_all <- learner$param_set
  res1 <- sapply(1:length(paras_to_tune), function(q){
    # q = 1
    tdf <- as.data.frame(as.data.table(ps_all))[which(paras_to_tune[q] == learner$param_set$ids()), ]
    if(tdf$class == "ParamInt"){
      if(tdf$upper == "Inf"){
        if(tdf$lower == "-Inf"){
          rts <- str_c(tdf$id, " = p_int(lower = ", -1e4, ", upper = ", 1e4, ")")  
        }else{
          rts <- str_c(tdf$id, " = p_int(lower = ", tdf$lower, ", upper = ", 1e4, ")")
        }
      }else{
        if(tdf$lower == "-Inf"){
          rts <- str_c(tdf$id, " = p_int(lower = ", -1e4, ", upper = ", tdf$upper, ")")
        }else{
          rts <- str_c(tdf$id, " = p_int(lower = ", tdf$lower, ", upper = ", tdf$upper, ")")
        }
      }
    }else if(tdf$class == "ParamDbl"){
      if(tdf$upper == "Inf"){
        if(tdf$lower == "-Inf"){
          rts <- str_c(tdf$id, " = p_dbl(lower = ", -1e4, ", upper = ", 1e4, ")")  
        }else{
          rts <- str_c(tdf$id, " = p_dbl(lower = ", tdf$lower, ", upper = ", 1e4, ")")
        }
      }else{
        if(tdf$lower == "-Inf"){
          rts <- str_c(tdf$id, " = p_dbl(lower = ", -1e4, ", upper = ", tdf$upper, ")")
        }else{
          rts <- str_c(tdf$id, " = p_dbl(lower = ", tdf$lower, ", upper = ", tdf$upper, ")")
        }
      }
    }else{
      rts <- ''
    }
    return(rts)
  })
  res1 <- res1[res1 != '']
  for (tid in 1:length(res1)) {
    if(tid != length(res1)){
      res1[tid] <- str_c(res1[tid], ',', collapse = '')  
    }
  }
  res1 <- c("search_space = ps(",res1, ")")
  res1 <- str_c(res1, collapse = '')
  return(res1)
}
#==== Constants ===========
##==== not changable ====
st1 <- read_csv(file = 'rawdata/A4BX6train_and_test_fill1.csv')
abbs <- read_csv(file = 'rawdata/Abbs.csv')
st <- clean_up(df = st1, abbs = abbs)
features <- st[, -ncol(st):(-ncol(st) + 1)] 
fea_new <- fea_slc_bycor(features, corr = 0.95)
means <- apply(fea_new, 2, mean)
sds <- apply(fea_new, 2, sd)
fea_new <- scale(fea_new, center = means, scale = sds)
bandgapd <- cbind(fea_new, bandgap = st$lumo - st$homo) 
bandgapd <- as.data.frame(bandgapd)
bandgapd$bandgap <- as.factor(as.numeric(bandgapd$bandgap <= 3.4))
NN <- nrow(bandgapd)
##==== changable ====
patition_ratio <- 0.7
data_for_ML <- bandgapd
#==== Main work flow ===============
## ==== SVM ====
### ==== temp constant ====
learner_to_use <- "classif.svm"
train_plot_title <- str_c(str_replace(learner_to_use, 
                                      pattern = 'classif.',
                                      replacement = ''),
                          ":train set", collapse = '')
test_plot_title <- str_c(str_replace(learner_to_use, 
                                     pattern = 'classif.', 
                                     replacement = ''),
                         ":test set", collapse = '')
### ==== workflow ====
task <- TaskClassif$new(id = "data_for_ML", 
                        backend = data_for_ML, 
                        target = tail(colnames(data_for_ML), 1))
splits <- partition(task, ratio = patition_ratio)
learner <- lrn(learner_to_use,
               predict_type = 'prob',
               type = "C-classification")
learner$param_set$set_values(kernel = "radial")
learner$param_set$set_values(cost = 1.53)
set_threads(learner, n = 10)
learner$train(task = task, row_ids = splits$train)
predictions_train <- learner$predict(task, row_ids = splits$train)
predictions_test <- learner$predict(task, row_ids = splits$test)
measures = msrs(c("classif.acc", "classif.auc"))
measures$score(predictions_test)
predictions_test$confusion
svm_train_plot <- predplot_classif(task = task, 
                                  learner = learner,
                                  prediction = predictions_train,
                                  type = 'roc',
                                  title = "SVM: train set",
                                  subtitle = str_c("AUC: ",
                                                   round(measures$classif.auc$score(predictions_train), 
                                                         digits = 3),
                                                   ", ACC: ", round(measures$classif.acc$score(predictions_train),
                                                                    digits = 3)))
ggsave(filename = "JAP/results/ML/SVM_train_ROC.jpeg", 
       plot = svm_train_plot, 
       width = unit(3, 'cm'), 
       height = unit(3, 'cm'))
svm_test_plot <- predplot_classif(task = task, 
                                  learner = learner,
                                  prediction = predictions_test,
                                  type = 'roc',
                                  title = "SVM: test set",
                                  subtitle = str_c("AUC: ",
                                                   round(measures$classif.auc$score(predictions_test), 
                                                         digits = 3),
                                                   ", ACC: ", round(measures$classif.acc$score(predictions_test),
                                                                    digits = 3)))
ggsave(filename = "JAP/results/ML/SVM_test_ROC.jpeg", 
       plot = svm_test_plot, 
       width = unit(3, 'cm'), 
       height = unit(3, 'cm'))

svm_train_plot <- predplot_classif(task = task, 
                                   learner = learner,
                                   prediction = predictions_train,
                                   type = 'prc',
                                   title = "SVM: train set",
                                   subtitle = str_c("AUC: ",
                                                    round(measures$classif.auc$score(predictions_train), 
                                                          digits = 3),
                                                    ", ACC: ", round(measures$classif.acc$score(predictions_train),
                                                                     digits = 3)))
ggsave(filename = "JAP/results/ML/SVM_train_PRC.jpeg", 
       plot = svm_train_plot, 
       width = unit(3, 'cm'), 
       height = unit(3, 'cm'))
svm_test_plot <- predplot_classif(task = task, 
                                  learner = learner,
                                  prediction = predictions_test,
                                  type = 'prc',
                                  title = "SVM: test set",
                                  subtitle = str_c("AUC: ",
                                                   round(measures$classif.auc$score(predictions_test), 
                                                         digits = 3),
                                                   ", ACC: ", round(measures$classif.acc$score(predictions_test),
                                                                    digits = 3)))
ggsave(filename = "JAP/results/ML/SVM_test_PRC.jpeg", 
       plot = svm_test_plot, 
       width = unit(3, 'cm'), 
       height = unit(3, 'cm'))

## ==== xgboost ====
### ==== temp constant ====
learner_to_use <- "classif.xgboost"
train_plot_title <- str_c(str_replace(learner_to_use, 
                                      pattern = 'regr.',
                                      replacement = ''),
                          ":train set", collapse = '')
test_plot_title <- str_c(str_replace(learner_to_use, 
                                     pattern = 'regr.', 
                                     replacement = ''),
                         ":test set", collapse = '')
### ==== workflow ====
task <- TaskClassif$new(id = "data_for_ML", 
                        backend = data_for_ML, 
                        target = tail(colnames(data_for_ML), 1))
splits <- partition(task, ratio = patition_ratio)
learner <- lrn(learner_to_use,
               predict_type = 'prob',
               nrounds = 6,
               booster = 'dart',
               subsample = 0.4,
               max_depth = 5,
               min_child_weight = 0,
               colsample_bytree = 1,
               colsample_bylevel = 0.01,
               rate_drop = 0,
               skip_drop = 0)
set_threads(learner, n = 10)
learner$train(task = task, row_ids = splits$train)
predictions_train <- learner$predict(task, row_ids = splits$train)
predictions_test <- learner$predict(task, row_ids = splits$test)
measure = msr("classif.acc")
measure$score(predictions_test)
measure = msr("classif.auc")
measure$score(predictions_test)
predictions_test$confusion

xgboost_train_plot <- predplot_classif(task = task, 
                                   learner = learner,
                                   prediction = predictions_train,
                                   type = 'roc',
                                   title = "xgboost: train set",
                                   subtitle = str_c("AUC: ",
                                                    round(measures$classif.auc$score(predictions_train), 
                                                          digits = 3),
                                                    ", ACC: ", round(measures$classif.acc$score(predictions_train),
                                                                     digits = 3)))
ggsave(filename = "JAP/results/ML/xgboost_train_ROC.jpeg", 
       plot = xgboost_train_plot, 
       width = unit(3, 'cm'), 
       height = unit(3, 'cm'))
xgboost_test_plot <- 
  predplot_classif(task = task, 
                   learner = learner,
                   prediction = predictions_test,
                   type = 'roc',
                   title = "xgboost: test set",
                   subtitle = str_c("AUC: ",
                                    round(measures$classif.auc$score(predictions_test), 
                                          digits = 3),
                                    ", ACC: ", 
                                    round(measures$classif.acc$score(predictions_test),
                                          digits = 3)))
ggsave(filename = "JAP/results/ML/xgboost_test_ROC.jpeg", 
       plot = xgboost_test_plot, 
       width = unit(3, 'cm'), 
       height = unit(3, 'cm'))

xgboost_train_plot <- predplot_classif(task = task, 
                                   learner = learner,
                                   prediction = predictions_train,
                                   type = 'prc',
                                   title = "SVM: train set",
                                   subtitle = str_c("AUC: ",
                                                    round(measures$classif.auc$score(predictions_train), 
                                                          digits = 3),
                                                    ", ACC: ", round(measures$classif.acc$score(predictions_train),
                                                                     digits = 3)))
ggsave(filename = "JAP/results/ML/xgboost_train_PRC.jpeg", 
       plot = xgboost_train_plot, 
       width = unit(3, 'cm'), 
       height = unit(3, 'cm'))
xgboost_test_plot <- predplot_classif(task = task, 
                                  learner = learner,
                                  prediction = predictions_test,
                                  type = 'prc',
                                  title = "xgboost: test set",
                                  subtitle = str_c("AUC: ",
                                                   round(measures$classif.auc$score(predictions_test), 
                                                         digits = 3),
                                                   ", ACC: ", round(measures$classif.acc$score(predictions_test),
                                                                    digits = 3)))
ggsave(filename = "JAP/results/ML/xgboost_test_PRC.jpeg", 
       plot = xgboost_test_plot, 
       width = unit(3, 'cm'), 
       height = unit(3, 'cm'))

##==== prediction ====







