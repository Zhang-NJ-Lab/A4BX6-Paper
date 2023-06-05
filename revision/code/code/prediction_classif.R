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
           -contains('Group'), -contains('Year'), -contains('NumberofValence'),
           -contains('name'))
  df <- df[!is.na(df$homo), ]
  colnames(df) <- abbs$Abbreviation
  df <- as.data.frame(df)
  st <- df[, -1]
  for(s in 1:ncol(st)){
    st[, s] <- as.numeric(st[, s])
  }
  return(st)
}
#==== Constants ===========
##==== not changable ====
st1 <- read_csv(file = 'rawdata/A4BX6train_and_test_fill1.csv')
abbs <- read_csv(file = 'rawdata/Abbs.csv')
st1$name <- str_c(st1$A1_Symbol,
                  st1$A2_Symbol,
                  st1$A3_Symbol,
                  st1$A4_Symbol,
                  st1$B_Symbol,
                  st1$X1_Symbol,
                  st1$X2_Symbol,
                  st1$X3_Symbol,
                  st1$X4_Symbol,
                  st1$X5_Symbol,
                  st1$X6_Symbol)
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
### ==== workflow ====
task <- TaskClassif$new(id = "data_for_ML", 
                        backend = data_for_ML, 
                        target = tail(colnames(data_for_ML), 1))
splits <- partition(task, ratio = patition_ratio)
train_indcies_tibble <- read_csv(file = 'JAP/results/train_index.csv')
test_indcies_tibble <- read_csv(file = 'JAP/results/test_index.csv')
splits$train <- train_indcies_tibble$train
splits$test <- test_indcies_tibble$test
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
measures$classif.acc$score(predictions_test)
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
# ggsave(filename = "JAP/results/ML/SVM_train_ROC.jpeg", 
#        plot = svm_train_plot, 
#        width = unit(3, 'cm'), 
#        height = unit(3, 'cm'))
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
# ggsave(filename = "JAP/results/ML/SVM_test_ROC.jpeg", 
#        plot = svm_test_plot, 
#        width = unit(3, 'cm'), 
#        height = unit(3, 'cm'))

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
# ggsave(filename = "JAP/results/ML/SVM_train_PRC.jpeg", 
#        plot = svm_train_plot, 
#        width = unit(3, 'cm'), 
#        height = unit(3, 'cm'))
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
# ggsave(filename = "JAP/results/ML/SVM_test_PRC.jpeg", 
#        plot = svm_test_plot, 
#        width = unit(3, 'cm'), 
#        height = unit(3, 'cm'))

### ==== prediciton on newdata ====
pred_set <- read_csv(file = 'MS文件生成/Results/prediction_cleaned.csv')
pred_set_dirty <- read_csv(file = 'MS文件生成/Results/prediction.csv')
pred_set <- pred_set[, colnames(data_for_ML)[1:(length(colnames(data_for_ML)) - 1)]]
pred_set <- as.data.frame(pred_set)
dim(pred_set)
pred_set <- scale(pred_set, center = means, scale = sds)
pred_set <- as.data.frame(pred_set)
learner$predict_type <- 'response'
pred_set_response <- learner$predict_newdata(newdata = pred_set)
pred_set$bandgap <- pred_set_response$response
pred_set$name <- pred_set_dirty$name
readable_pred_set <- data.frame(name = pred_set$name,
                                solvent = pred_set$solvent,
                                bandgap = pred_set$bandgap)
tmp <- data.frame(solvent = scale(st1$solvent, center = means['solvent'], scale = sds['solvent']))
st1$solvent <- tmp$solvent
name_extract <- readable_pred_set$name %in% st1$name
res2 <- mclapply(1:nrow(readable_pred_set[name_extract, ]), function(q){
  df <- readable_pred_set[name_extract, ][q, 'solvent'] == st1[st1$name == readable_pred_set[name_extract, ][q, 'name'],'solvent']
  return(df[1, 1])
}, mc.cores = 20)
res2 <- unlist(res2)
readable_pred_set1 <- readable_pred_set[(!res2), ]
tmp <- data.frame(X1 = scale(sol$X1, 
                                  center = means['solvent'], 
                                  scale = sds['solvent']))
sol$X1 <- tmp$X1
for (s in 1:nrow(readable_pred_set1)) {
  readable_pred_set1$solvent[s] <- sol$X2[readable_pred_set1[s, 'solvent'] == sol$X1]
}
# write_csv(x = pred_set, file = 'JAP/results/preditcion_raw.csv')
write_csv(x = readable_pred_set1, file = 'JAP/results/readable_pred_set.csv')
write_csv(x = readable_pred_set[readable_pred_set1$bandgap !=0, ], file = 'JAP/results/readable_pred_set_positive.csv')


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







