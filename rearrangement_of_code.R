rm(list = ls()) 
setwd("~/Desktop/ZhengPan/workspace/NUIST/John/胡文广/A4BX6/") 
source('code/dependencies.R')
#= data import ============
st1 <- read_csv(file = 'rawdata/A4BX6train_and_test_fill1.csv') 
st1 <- st1 %>%
  select(-contains('Symbol'), -contains('Phase'), -contains('Radioactive'), 
         -contains('Natural'), -contains('Metal'), -contains('Nonmetal'), 
         -contains('Metalloid'), -contains('Type'),-contains('Discoverer'),
         -contains('Group'), -contains('Year'), -contains('NumberofValence'))
st1 <- st1[!is.na(st1$homo), ]
abbs <- read_csv(file = 'rawdata/Abbs.csv')
colnames(st1) <- abbs$Abbreviation
st1 <- as.data.frame(st1)
st <- st1[, -1]
for(s in 1:ncol(st)){
  st[, s] <- as.numeric(st[, s])
}
#### main #### 
NN <- nrow(st1)
# feature engineering ## Pearson heat map colnames(st)
dim(st)
features <- st[, -ncol(st):(-ncol(st) + 1)] 
colnames(features)
rawFeaPic <- pearson_heat(features, corm = cor(features)) + 
  theme(legend.position = "right",
        legend.title = element_text(size = 10), 
        legend.key.height = unit(2, "cm"),
        axis.text.x = element_text(size = 3.3,
                                   angle = 90, 
                                   hjust = .5, 
                                   vjust = .4),
        axis.text.y = element_text(size = 3.3,
                                   angle = 0)
  ) + 
  scale_x_discrete(labels = fancy_axis) + 
  scale_y_discrete(labels = fancy_axis)
# ggsave(filename = "results/feature_engineering/Raw_pearson_heat.jpeg", plot = rawFeaPic, width = 17.24, height = 13.41, units = 'cm')
fea_new <- fea_slc_bycor(features, corr = 0.95)
homod <- cbind(fea_new, homo = st$homo)
lumod <- cbind(fea_new, lumo = st$lumo)
bandgapd <- cbind(fea_new, bandgap = st$lumo - st$homo)
fea_new_names <- colnames(fea_new)
newFeaPic <- pearson_heat(fea_new) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 20), 
        legend.key.height = unit(3, "cm"),
        # text = element_text(size = 30), 
        axis.text.x = element_text(size = 18,
                                   angle = 270 + 45, hjust = .2,
                                   vjust = .3),
        axis.text.y = element_text(size = 18, angle = 0)) +
  scale_x_discrete(labels = fancy_axis) +
  scale_y_discrete(labels = fancy_axis)
# ggsave(filename = "results/feature_engineering/ Feature_selected_pearson_heat.jpeg", plot = newFeaPic)
## correlation
cor_for_fea_sel_df_pic_homo <- cwg(fea_new_andC = homod, lim = c(-0.6, 0.6), hlab = 'homo')
cor_for_fea_sel_df_pic_lumo <- cwg(fea_new_andC = lumod, lim = c(-0.6, 0.5), hlab = 'lumo')
cor_for_fea_sel_df_pic_bandgap <- cwg(fea_new_andC = bandgapd, lim = c(-0.6, 0.3), hlab = 'bandgap')
# ggsave(filename = "results/feature_engineering/feature_ranking_homo.jpeg", plot = cor_for_fea_sel_df_pic_homo, width = unit(7,"in"), height = unit(10,"in"))
# ggsave(filename = "results/feature_engineering/feature_ranking_lumo.jpeg", plot = cor_for_fea_sel_df_pic_lumo, width = unit(7,"in"), height = unit(10,"in"))
# ggsave(filename = "results/feature_engineering/feature_ranking_bandgap.jpeg", plot = cor_for_fea_sel_df_pic_bandgap, width = unit(7,"in"), height = unit(10,"in"))

# GE =========================


#= symbolic regression ======================
# test1
ipd <- fea_new
opd <- homod$homo
da <- cbind(ipd,opd)
fn <- names(ipd)
fn <- colnames(fea_new)
tn <- names(opd)
da <- as.data.frame(da)
da1 <- as.data.frame(matrix(0, nr = nrow(da), nc = ncol(da)))
for (s in 1:nrow(da)) {
  da1[s, ] <- unlist(da[s, ])
}
da <- da1
da <- as.data.frame(scale(da))
colnames(da) <- c(paste("feature",1:(ncol(da)-1),sep = "_"),"tgt")
colnames(da1) <- c(paste("feature",1:(ncol(da1)-1),sep = "_"),"tgt")
attach(da)
ruleDef <- list(expr = grule(op(expr, expr), func(expr), var),
                func = grule(log, sqrt),
                op = grule(`+`, `-`, `*`, `/`, `^`),
                op = grule(`+`, `-`, `*`, `/`),
                var = grule(feature_1^n, feature_2^n, feature_3^n, feature_4^n, feature_5^n, feature_6^n, feature_7^n, feature_8^n, feature_9^n, feature_10^n, feature_11^n, feature_12^n, feature_13^n),
                n = grule(1, 2, 3)
)
grammarDef <- CreateGrammar(ruleDef)
print(grammarDef)

SymRegFitFunc <- function(expr) { 
  result <- eval(expr)
  if (any(is.nan(result)) | any(is.na(result))){
    return(Inf)
  }
  if(is.nan((cor(cbind(unlist(result), unlist(opd)))[1,2]))){
    return(Inf)
  }
  if(is.na((cor(cbind(unlist(result), unlist(opd)))[1,2]))){
    return(Inf)
  }
  if((cor(cbind(unlist(result), unlist(opd)))[1,2]) <= 0){
    return(Inf)
  }
  return(1/abs(cor(result, opd)))
}
library(parallel)
ge <- GrammaticalEvolution(grammarDef = grammarDef, 
                           evalFunc = SymRegFitFunc, 
                           terminationCost = 1.001,
                           monitorFunc = print,
                           optimizer = "ga",
                           max.depth = 5,
                           seqLen = 50,
                           iterations = 1e5,
                           popSize = 1e3,
                           mutationChance = 0.3,
                           plapply = mclapply
)

# ACM-AMI_1 
## without scale
best.expression <- feature_3^2 * feature_11^1 * abs(feature_11^3 - feature_9^3)
best.expression <- abs(sqrt(feature_11^3) - feature_9^2) 
best.expression <- abs(abs(feature_11^3/abs(abs(abs(feature_9^3)) * feature_12^3))) + feature_9^3 
best.expression <- (feature_11^3 - feature_9^3) * (sqrt(feature_8^1) - feature_9^2) 
best.expression <- abs(feature_10^3/(log(log(feature_9^3) * feature_11^1) - sqrt(feature_10^1) - abs(feature_6^1))) 
## with scale
best.expression <- feature_7^1 - feature_12^1 - feature_4^1 - (feature_6^1 + log(feature_11^2 + feature_4^2 - feature_4^2) - (sqrt(feature_2^2) + feature_9^1))

cor(as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))))[1,2]
as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))) %>%
  ggplot(aes(x = x)) +
  geom_point(aes(y = opd), color = 'red') +
  geom_point(aes(y = GE), color = 'blue') +
  theme_bw() +
  labs(title = "(cos(sin(0.8787879 * feature_5^1) + cos(0.959596 * feature_12^1.210526)) + 0.5151515 * \nfeature_5^1)/(0.3737374 * feature_12^1.842105) * (0.1414141 * feature_9^3.315789) ")

as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))) %>%
  ggplot(aes(x = GE, y = opd)) +
  geom_point(color = 'red') +
  geom_smooth(method = 'lm',formula = y ~ x - 1) +
  theme_bw()
detach(da)
ggsave(filename = '~/desktop/ge/test2.jpeg')

da <- as.data.frame(da)
attach(da)

detach(da)
feature_1

expr_without_scale <- list(
  expr(feature_3^2 * feature_11^1 * abs(feature_11^3 - feature_9^3)),
  expr(abs(sqrt(feature_11^3) - feature_9^2)) ,
  expr(abs(abs(feature_11^3/abs(abs(abs(feature_9^3)) * feature_12^3))) + feature_9^3 ),
  expr((feature_11^3 - feature_9^3) * (sqrt(feature_8^1) - feature_9^2) ),
  expr(abs(feature_10^3/(log(log(feature_9^3) * feature_11^1) - sqrt(feature_10^1) - abs(feature_6^1)))) 
)
expr_with_scale <- list(
  expr(feature_9^1 * (feature_3^2 - feature_11^3)),
  expr((feature_3^1 - feature_11^1) * feature_9^1),
  expr(feature_3^3/(feature_11^2 + feature_8^2)),
  expr((feature_9^3 - sqrt(sqrt(feature_4^2))) * (feature_3^1 - feature_11^1)),
  expr(feature_3^2 - feature_12^1),
  expr(feature_3^2 - feature_9^1 * feature_11^1 ),
  expr(feature_9^2 - (feature_11^1 - feature_6^2) * feature_9^1 ),
  expr(feature_3^2/((feature_11^2 + feature_4^2/(feature_9^2/feature_10^2)) * feature_5^2) ),
  expr(feature_4^2 * (feature_3^2 + feature_11^3/(feature_11^2 + feature_10^1)) ),
  expr(log(feature_3^2) * (feature_3^3 - feature_11^3 - feature_8^3) ),
  expr(feature_9^1 * (feature_9^2 - (feature_11^1 + sqrt(feature_12^2))) ),
  expr(feature_1^1/(feature_7^1 * (feature_11^2 - feature_4^1)) * feature_3^2 )
)



attach(da1)
corWithoutScale <- lapply(1:length(expr_without_scale), function(q){
  best.expression <- expr_without_scale[q]
  cor(as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))))[1,2]
})
detach(da1)
corWithoutScale

df1 <- data.frame(formula = c(sapply(1:length(corWithoutScale), function(q){
  str_c('F', q + 12)
}), rownames(cor_for_fea_sel_df)), correlation = c(unlist(corWithoutScale), unlist(cor_for_fea_sel_df))) %>%
  arrange((correlation))
ggplot(df1) +
  geom_col(aes(y = correlation, x = factor(formula, levels = (formula))), fill = "yellow") +
  theme_bw() +
  coord_flip() +
  # scale_y_continuous(limits = c(-0.2,0.4)) +
  xlab("Features") +
  ylab("Pearson Coefficient") +
  theme(axis.text.y = element_text(size = 15,
                                   angle = 0,
                                   hjust = 0.8,
                                   vjust = 0.0),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.margin = margin(10, 10, 10, 10)) +
  # scale_x_continuous(labels = fancy_axis) +
  scale_x_discrete(labels = fancy_axis) +
  annotate('text',
           x = df1$formula,
           y = df1$correlation, 
           label = round(df1$correlation, 2))
ggsave(filename = 'ACS AMI/pz1/NotScale.jpeg', width = unit(7,"in"), height = unit(10,"in"))


attach(da)
corWithScale <- lapply(1:length(expr_with_scale), function(q){
  best.expression <- expr_with_scale[q]
  cor(as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))))[1,2]
})
detach(da)

df2 <- data.frame(formula = c(sapply(1:length(corWithScale), function(q){
  str_c('F', q)
}), rownames(cor_for_fea_sel_df)), correlation = c(unlist(corWithScale), unlist(cor_for_fea_sel_df))) %>%
  arrange((correlation))
ggplot(df2) +
  geom_col(aes(y = correlation, x = factor(formula, levels = (formula))), fill = "yellow") +
  theme_bw() +
  coord_flip() +
  # scale_y_continuous(limits = c(-0.2,0.4)) +
  xlab("Features") +
  ylab("Pearson Coefficient") +
  theme(axis.text.y = element_text(size = 15,
                                   angle = 0,
                                   hjust = 0.8,
                                   vjust = 0.0),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.margin = margin(10, 10, 10, 10)) +
  # scale_x_continuous(labels = fancy_axis) +
  scale_x_discrete(labels = fancy_axis) +
  annotate('text',
           x = df2$formula,
           y = df2$correlation, 
           label = round(df2$correlation, 2))
ggsave(filename = 'ACS AMI/pz1/Scale.jpeg', width = unit(7,"in"), height = unit(10,"in"))

## change split ratio of KM ============
task <- TaskRegr$new(id = "task", backend = tsk_df[1:NN, ], target = "current")
split_ratio <- seq(0.1, 1, 0.1)
cl <- mlr_learners$keys()[35:44]
qq = 5
learner <- mlr_learners$get(cl[qq])
res1 <- lapply(16:(NN) - 1, function(q){
  # q = 16
  train_set <- sample(1:NN, size = q)
  test_set<- setdiff(1:NN, train_set)
  learner$train(task, row_ids = train_set)
  prediction <- learner$predict(task, row_ids = test_set)
  measure <- msr("regr.mse")
  mse <- format(prediction$score(measure), digits = 5)
  measure <- msr("regr.rmse")
  rmse <- format(prediction$score(measure), digits = 2)
  preds <- prediction$response
  actual <- prediction$truth
  pr <- round(cor(cbind(preds, actual))[1,2],2)
  rss <- sum((preds - actual) ^ 2)  ## residual sum of squares
  tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
  rsq <- round(1 - rss/tss, 2)
  # test_set_pic_km <- predplot(prediction, title = 'Kriging model: test set')
  return(list(pr, rmse, rsq))
})
res2 <- do.call(rbind, res1)
res2 <- as.data.frame(res2)
write.csv(res2, file = 'temp1.csv')
res2 <- read_csv(file = 'temp1.csv')
colnames(res2) <- c('i','pr', 'rmse', 'rsq')

ggplot(data = res2) +
  geom_point(aes(x = i, y = rmse), size = 0.5) +
  geom_smooth(aes(x = i, y = rmse),se = F, color = 'red') +
  theme_bw() +
  xlab('Training set') +
  ylab('RMSE') +
  scale_y_continuous(labels = fancy_scientific)
ggsave(filename = "ACS AMI/pics/ratio.jpeg", width = unit(3,"cm"), height = unit(3,"cm"))


## === pearson head map with F ==========

attach(da1)
corWithoutScale <- lapply(1:length(expr_without_scale), function(q){
  # q = 1
  best.expression <- expr_without_scale[q][[1]]
  return(eval(best.expression))
  # cor(as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))))[1,2]
})
detach(da1)
F2 <- do.call(cbind, corWithoutScale)
colnames(F2) <- sapply(1:length(corWithoutScale), function(q){
  str_c('F', q + 12)
})


attach(da)
corWithScale <- lapply(1:length(expr_with_scale), function(q){
  best.expression <- expr_with_scale[q][[1]]
  return(eval(best.expression))
})
detach(da)
F1 <- do.call(cbind, corWithScale)
colnames(F1) <- sapply(1:length(corWithScale), function(q){
  str_c('F', q)
})

pearson_heat(cbind(fea_new, F1, F2)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        legend.key.height = unit(2, "cm"),
        # text = element_text(size = 30),
        axis.text.x = element_text(
          size = 7,
          angle = 45 + 90 + 90 + 90,
          hjust = .1,
          vjust = .9),
        axis.text.y = element_text(
          size = 7,
          angle = 0)) +
  scale_x_discrete(labels = fancy_axis) +
  scale_y_discrete(labels = fancy_axis)
ggsave(filename = "ACS AMI/pics/newp.jpeg", width = unit(5,"cm"), height = unit(5,"cm"))








#### test ####=========================
#### visualization of prediction ####
st_t <- read_csv(file = "results/predictions/KMpredions.csv")
st_t <- as.data.frame(st_t)
st_t_2 <- read_csv(file = "rawdata/光电流数据yingwenred.csv")
st_t_2 <- as.data.frame(st_t_2)
colnames(st_t_2) <- colnames(st_t)
st_t <- rbind(st_t_2, st_t)
st_t$lab <- "a"
st_t$lab[(NN + 1): nrow(st_t)] <- "prediction set"
km_test_set <- c(unlist(as.data.frame(read_csv(file = "results/temp_data/ test_set.csv"))))
names(km_test_set) <- NULL
km_train_set <- unlist(as.data.frame(read_csv(file = "results/temp_data/ train_set.csv")))
names(km_train_set) <- NULL
st_t$lab[km_test_set] <- "test set"
st_t$lab[km_train_set] <- "train set"
st_t_1 <- rbind(st_t[st_t$lab == "prediction set", ], st_t[st_t$lab == "train set", ], st_t[st_t$lab == "test set", ])

colnames(st_t_1)
ggplot(data = st_t_1[st_t_1$Current.A. <= 5e-6 &st_t_1$Current.A. >= 4.75e-6 , ], mapping = aes(y = `Current.A.`))+
  geom_path(aes(x = `Time.s.`, group = factor(X), color = factor(lab), fill = factor(lab)), alpha = 0.5) +
  theme_bw() + scale_color_manual(
    name = NULL,
    breaks = c("prediction set", "train set", 'test set'), values = c(
      'prediction set' = '#00FFFF', 'train set' = '#FF0000',
      'test set' = '#000066'
    ) )
st_t_1$lab <- factor(st_t_1$lab, levels = sort(unique(st_t_1$lab))) p_d_1 <- ggplot(data = st_t_1, mapping = aes(y = `Current.A.`))+
  geom_point(aes(x = `Time.s.`,color = factor(lab)), alpha = 0.4, size = .4) + theme_bw() +
  scale_color_manual(
    name = NULL,
    breaks = c("prediction set", "train set", 'test set'), values = c(
      'prediction set' = '#00FFFF', 'train set' = '#FF0000',
      'test set' = '#000066'
    ) )+
  theme(
    legend.position = c(.92, .92),
    legend.justification = c(1,1),
    legend.background = element_rect(colour = 'blue', fill = '#FFFFCC'), legend.key = element_rect(colour = NULL, fill = '#FFFFCC'),
    text = element_text(size = 15)
  )+
  guides(color = guide_legend(title=NULL)) + ylab('Photocurrent (A)') +
  xlab('Time (s)')
ggsave(p_d_1, filename = "results/pred_visualization/p_d_1.jpeg", width = 5.32, height = 5.08)

p_d_2 <- ggplot(data = st_t_1, mapping = aes(y = `Current.A.`))+ geom_point(aes(x = `Time.s.`,color = factor(X)), alpha = 0.3, size = .4) + theme_bw() +
  theme(
    legend.position = c(.92, .92),
    legend.justification = c(1,1),
    legend.background = element_rect(colour = 'blue', fill = '#FFFFCC'), legend.key.size = unit(.3, 'cm'),
    legend.key = element_rect(colour = NULL, fill = '#FFFFCC'), axis.title = element_text(size = 15),
    axis.text = element_text(size = 15)
  )+
  guides(color = guide_legend(title=NULL)) + ylab('Photocurrent (A)') +
  xlab('Time (s)')
ggsave(p_d_2, filename = "results/pred_visualization/p_d_2.jpeg", width = 5.32, height = 5.08)
#### plot for 12 features ####
plot_of_prediation <- function(f_df, f_xn, f_yn){ # f_df = st_t_1
  # f_xn = 1
  # f_yn = ncol(st_t_1) - 1
  #t=1 res_p <-
  ggplot(data = f_df, mapping = aes(y = f_df[, f_yn])) +
    # geom_point(aes(x = jitter(f_df[, fea_new_names[f_xn]], amount = .005), color = factor(lab)), alpha = 0.4, size = .5) +
    geom_point(aes(x = f_df[, fea_new_names[f_xn]], color = factor(lab)), alpha = 0.4, size = .5) +
    theme_bw() + scale_color_manual(
      name = NULL,
      breaks = c("prediction set", "train set", 'test set'), values = c(
        'prediction set' = '#00FFFF', 'train set' = '#FF0000',
        'test set' = '#000066'
      )
      
    )+
    xlab(fea_new_names[f_xn]) +
    ylab('Current (A)') +
    guides(color = guide_legend(title=NULL)) + theme(
      legend.position = c(.92, .92),
      legend.justification = c(1,1),
      legend.background = element_rect(colour = 'blue', fill = '#FFFFCC'), legend.key.size = unit(.1, 'cm'),
      legend.key = element_rect(colour = NULL, fill = '#FFFFCC'), axis.title = element_text(size = 15)
    )+
    scale_y_continuous(labels = fancy_scientific) + xlab(fancy_axis(fea_new_names[f_xn]))
  # print(res_p)
  # ggsave(res_p, filename = paste0("results/pred_visualization/real", fea_new_names[t], ".jpeg"), width = 7.6, height = 7.9, units = 'cm')
  return(res_p) }
colnames(st_t_1)[1:(ncol(st_t_1) - 1)] <- colnames(st1)
for (t in 1:length(fea_new_names)) { t=8
print(paste0("starting job: ", t, " total job amount: ", length(fea_new_names), ", ", round(t/length(fea_new_names), 3)*100, " %"))
t_pic <- plot_of_prediation(st_t_1, t, ncol(st_t_1) - 1) temp_pic <- t_pic +
  theme(legend.position = c(.5, .92))
ggsave(temp_pic, filename = paste0("results/pred_visualization/x=",
                                   fea_new_names[t], ".jpeg"), width = 7.9, height = 7.6, units = 'cm') ggsave(t_pic, filename = paste0("results/pred_visualization/x=", fea_new_names[t], ".jpeg"), width = 7.9, height = 7.6, units = 'cm')
}
fea_new_names[1]
geom_point(aes(x = `Time.s.`,color = factor(X)), alpha = 0.3, size = .4) +
  theme(
    legend.position = c(.92, .92),
    
    legend.justification = c(1,1),
    legend.background = element_rect(colour = 'blue', fill = '#FFFFCC'), legend.key.size = unit(.1, 'cm'),
    legend.key = element_rect(colour = NULL, fill = '#FFFFCC'), axis.title = element_text(size = 15)
  )+
  guides(color = guide_legend(title=NULL)) + ylab('Photocurrent (A)') +
  xlab('Time (s)')
ggsave(p_d_2, filename = "results/pred_visualization/p_d_2.jpeg", width = 5.32, height = 5.08)
# === AMI reply letter ===========================
#= Model Interpretation======= #== IML ========
library("iml")
x = penguins[which(names(penguins) != "species")]
model = Predictor$new(learner, data = x, y = penguins$species) num_features = c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g", "year")
effect = FeatureEffects$new(model)
plot(effect, features = num_features)
x = tsk_df[, which(colnames(tsk_df) != 'photocurrent')]
model = Predictor$new(learner, data = x, y = tsk_df$photocurrent) effect = FeatureEffects$new(model)
do.call(rbind, effect$results) %>%
  ggplot() +
  geom_line(mapping = aes(x = .borders, y = .value)) + facet_wrap(~.feature, ncol = 4, scales = "free_x") +
  
  theme_bw()
#= symbolic regression ======================
# test1
ipd <- fea_new
opd <- data.frame(opd = read.csv(file = "rawdata/光电流数据 yingwenred.csv", header = T)[, ncol(read.csv(file = "rawdata/光电流数据 yingwenred.csv", header = T))])
da <- cbind(ipd,opd)
fn <- names(ipd)
fn <- colnames(fea_new)
tn <- names(opd)
da <- as.data.frame(da)
da1 <- as.data.frame(matrix(0, nr = nrow(da), nc = ncol(da))) for (s in 1:nrow(da)) {
  da1[s, ] <- unlist(da[s, ]) }
da <- da1
da <- as.data.frame(scale(da))
colnames(da) <- c(paste("feature",1:(ncol(da)-1),sep = "_"),"tgt") colnames(da1) <- c(paste("feature",1:(ncol(da1)-1),sep = "_"),"tgt") attach(da)
ruleDef <- list(expr = grule(op(expr, expr), func(expr), var),
                func = grule(log, sqrt),
                op = grule(`+`, `-`, `*`, `/`, `^`),
                op = grule(`+`, `-`, `*`, `/`),
                var = grule(feature_1^n, feature_2^n, feature_3^n, feature_4^n,
                            feature_5^n, feature_6^n, feature_7^n, feature_8^n, feature_9^n, feature_10^n, feature_11^n, feature_12^n, feature_13^n),
                n = grule(1, 2, 3)
)
grammarDef <- CreateGrammar(ruleDef) print(grammarDef)
SymRegFitFunc <- function(expr) {
  result <- eval(expr)
  if (any(is.nan(result)) | any(is.na(result))){
    return(Inf) }
  if(is.nan((cor(cbind(unlist(result), unlist(opd)))[1,2]))){
    
    return(Inf) }
  if(is.na((cor(cbind(unlist(result), unlist(opd)))[1,2]))){ return(Inf)
  }
  if((cor(cbind(unlist(result), unlist(opd)))[1,2]) <= 0){
    return(Inf) }
  return(1/abs(cor(result, opd))) }
library(parallel)
ge <- GrammaticalEvolution(grammarDef = grammarDef,
)
evalFunc = SymRegFitFunc, terminationCost = 2.0001, monitorFunc = print, optimizer = "ga",
max.depth = 5, seqLen = 50, iterations = 1e5, popSize = 1e3, mutationChance = 0.3, plapply = mclapply
# ACM-AMI_1
## without scale
best.expression <- feature_3^2 * feature_11^1 * abs(feature_11^3 - feature_9^3)
best.expression <- abs(sqrt(feature_11^3) - feature_9^2)
best.expression <- abs(abs(feature_11^3/abs(abs(abs(feature_9^3)) * feature_12^3))) + feature_9^3
best.expression <- (feature_11^3 - feature_9^3) * (sqrt(feature_8^1) - feature_9^2)
best.expression <- abs(feature_10^3/(log(log(feature_9^3) * feature_11^1) - sqrt(feature_10^1) - abs(feature_6^1)))
## with scale
best.expression <- feature_9^1 * (feature_3^2 - feature_11^3) best.expression <- (feature_3^1 - feature_11^1) * feature_9^1 best.expression <- feature_3^3/(feature_11^2 + feature_8^2) best.expression <- (feature_9^3 - sqrt(sqrt(feature_4^2))) * (feature_3^1 - feature_11^1)
best.expression <- feature_3^2 - feature_12^1

feature_3^2 - feature_9^1 * feature_11^1
feature_9^2 - (feature_11^1 - feature_6^2) * feature_9^1 feature_3^2/((feature_11^2 + feature_4^2/(feature_9^2/feature_10^2)) * feature_5^2)
feature_4^2 * (feature_3^2 + feature_11^3/(feature_11^2 + feature_10^1)) feature_3^2 - feature_12
log(feature_3^2) * (feature_3^3 - feature_11^3 - feature_8^3)
feature_9^1 * (feature_9^2 - (feature_11^1 + sqrt(feature_12^2))) feature_1^1/(feature_7^1 * (feature_11^2 - feature_4^1)) * feature_3^2
cor(as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))))[1,2] as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))) %>%
  ggplot(aes(x = x)) +
  geom_point(aes(y = opd), color = 'red') +
  geom_point(aes(y = GE), color = 'blue') +
  theme_bw() +
  labs(title = "(cos(sin(0.8787879 * feature_5^1) + cos(0.959596 *
feature_12^1.210526)) + 0.5151515 * \nfeature_5^1)/(0.3737374 * feature_12^1.842105) * (0.1414141 * feature_9^3.315789) ")
as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))) %>%
  ggplot(aes(x = GE, y = opd)) +
  geom_point(color = 'red') + geom_smooth(method = 'lm',formula = y ~ x - 1) + theme_bw()
detach(da)
ggsave(filename = '~/desktop/ge/test2.jpeg')
da <- as.data.frame(da) attach(da)
detach(da) feature_1
expr_without_scale <- list(
  expr(feature_3^2 * feature_11^1 * abs(feature_11^3 - feature_9^3)), expr(abs(sqrt(feature_11^3) - feature_9^2)) , expr(abs(abs(feature_11^3/abs(abs(abs(feature_9^3)) * feature_12^3))) +
                                                                                                                           feature_9^3 ),
  
  expr((feature_11^3 - feature_9^3) * (sqrt(feature_8^1) - feature_9^2) ),
  expr(abs(feature_10^3/(log(log(feature_9^3) * feature_11^1) - sqrt(feature_10^1) - abs(feature_6^1))))
)
expr_with_scale <- list(
  expr(feature_9^1 * (feature_3^2 - feature_11^3)),
  expr((feature_3^1 - feature_11^1) * feature_9^1), expr(feature_3^3/(feature_11^2 + feature_8^2)),
  expr((feature_9^3 - sqrt(sqrt(feature_4^2))) * (feature_3^1 - feature_11^1)), expr(feature_3^2 - feature_12^1),
  expr(feature_3^2 - feature_9^1 * feature_11^1 ),
  expr(feature_9^2 - (feature_11^1 - feature_6^2) * feature_9^1 ), expr(feature_3^2/((feature_11^2 + feature_4^2/(feature_9^2/feature_10^2))
                                                                                     * feature_5^2) ),
  expr(feature_4^2 * (feature_3^2 + feature_11^3/(feature_11^2 +
                                                    feature_10^1)) ),
  expr(log(feature_3^2) * (feature_3^3 - feature_11^3 - feature_8^3) ), expr(feature_9^1 * (feature_9^2 - (feature_11^1 + sqrt(feature_12^2))) ), expr(feature_1^1/(feature_7^1 * (feature_11^2 - feature_4^1)) *
                                                                                                                                                         feature_3^2 ) )
attach(da1)
corWithoutScale <- lapply(1:length(expr_without_scale), function(q){
  best.expression <- expr_without_scale[q]
  cor(as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))))[1,2]
})
detach(da1) corWithoutScale
df1 <- data.frame(formula = c(sapply(1:length(corWithoutScale), function(q){ str_c('F', q + 12)
}), rownames(cor_for_fea_sel_df)), correlation = c(unlist(corWithoutScale), unlist(cor_for_fea_sel_df))) %>%
  arrange((correlation)) ggplot(df1) +
  geom_col(aes(y = correlation, x = factor(formula, levels = (formula))), fill = "yellow") +
  
  theme_bw() +
  coord_flip() +
  # scale_y_continuous(limits = c(-0.2,0.4)) + xlab("Features") +
  ylab("Pearson Coefficient") + theme(axis.text.y = element_text(size = 15,
                                                                 angle = 0, hjust = 0.8, vjust = 0.0),
                                      axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), plot.margin = margin(10, 10, 10, 10)) +
  # scale_x_continuous(labels = fancy_axis) +
  scale_x_discrete(labels = fancy_axis) + annotate('text',
                                                   x = df1$formula,
                                                   y = df1$correlation,
                                                   label = round(df1$correlation, 2))
ggsave(filename = 'ACS AMI/pz1/NotScale.jpeg', width = unit(7,"in"), height = unit(10,"in"))
attach(da)
corWithScale <- lapply(1:length(expr_with_scale), function(q){
  best.expression <- expr_with_scale[q]
  cor(as.data.frame(cbind(data.frame(opd = opd, GE = eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))))[1,2]
})
detach(da)
df2 <- data.frame(formula = c(sapply(1:length(corWithScale), function(q){ str_c('F', q)
}), rownames(cor_for_fea_sel_df)), correlation = c(unlist(corWithScale), unlist(cor_for_fea_sel_df))) %>%
  arrange((correlation)) ggplot(df2) +
  geom_col(aes(y = correlation, x = factor(formula, levels = (formula))), fill = "yellow") +
  theme_bw() + coord_flip() +
  
  # scale_y_continuous(limits = c(-0.2,0.4)) + xlab("Features") +
  ylab("Pearson Coefficient") + theme(axis.text.y = element_text(size = 15,
                                                                 angle = 0, hjust = 0.8, vjust = 0.0),
                                      axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), plot.margin = margin(10, 10, 10, 10)) +
  # scale_x_continuous(labels = fancy_axis) +
  scale_x_discrete(labels = fancy_axis) + annotate('text',
                                                   x = df2$formula,
                                                   y = df2$correlation,
                                                   label = round(df2$correlation, 2))
ggsave(filename = 'ACS AMI/pz1/Scale.jpeg', width = unit(7,"in"), height = unit(10,"in"))
## change split ratio of KM ============
task <- TaskRegr$new(id = "task", backend = tsk_df[1:NN, ], target = "current")
split_ratio <- seq(0.1, 1, 0.1)
cl <- mlr_learners$keys()[35:44]
qq = 5
learner <- mlr_learners$get(cl[qq])
res1 <- lapply(16:(NN) - 1, function(q){
  # q = 16
  train_set <- sample(1:NN, size = q)
  test_set<- setdiff(1:NN, train_set)
  learner$train(task, row_ids = train_set)
  prediction <- learner$predict(task, row_ids = test_set)
  measure <- msr("regr.mse")
  mse <- format(prediction$score(measure), digits = 5)
  measure <- msr("regr.rmse")
  rmse <- format(prediction$score(measure), digits = 2)
  preds <- prediction$response
  actual <- prediction$truth
  pr <- round(cor(cbind(preds, actual))[1,2],2)
  rss <- sum((preds - actual) ^ 2) ## residual sum of squares
  tss <- sum((actual - mean(actual)) ^ 2) ## total sum of squares
  
  rsq <- round(1 - rss/tss, 2)
  # test_set_pic_km <- predplot(prediction, title = 'Kriging model: test set') return(list(pr, rmse, rsq))
})
res2 <- do.call(rbind, res1)
res2 <- as.data.frame(res2) write.csv(res2, file = 'temp1.csv')
res2 <- read_csv(file = 'temp1.csv') colnames(res2) <- c('i','pr', 'rmse', 'rsq')
ggplot(data = res2) +
  geom_point(aes(x = i, y = rmse), size = 0.5) + geom_smooth(aes(x = i, y = rmse),se = F, color = 'red') + theme_bw() +
  xlab('Training set') +
  ylab('RMSE') +
  scale_y_continuous(labels = fancy_scientific)
ggsave(filename = "ACS AMI/pics/ratio.jpeg", width = unit(3,"cm"), height = unit(3,"cm"))
## === pearson head map with F ==========
attach(da1)
corWithoutScale <- lapply(1:length(expr_without_scale), function(q){
  #q=1
  best.expression <- expr_without_scale[q][[1]] return(eval(best.expression))
  # cor(as.data.frame(cbind(data.frame(opd = opd, GE =
  eval(best.expression)), x = 1:nrow(data.frame(opd = tgt, GE = eval(best.expression))))))[1,2]
})
detach(da1)
F2 <- do.call(cbind, corWithoutScale)
colnames(F2) <- sapply(1:length(corWithoutScale), function(q){
  str_c('F', q + 12) })
attach(da)
corWithScale <- lapply(1:length(expr_with_scale), function(q){
  best.expression <- expr_with_scale[q][[1]]
  
  return(eval(best.expression)) })
detach(da)
F1 <- do.call(cbind, corWithScale)
colnames(F1) <- sapply(1:length(corWithScale), function(q){
  str_c('F', q) })
pearson_heat(cbind(fea_new, F1, F2)) + theme(legend.position = "right",
                                             legend.text = element_text(size = 10), legend.key.height = unit(2, "cm"),
                                             # text = element_text(size = 30), axis.text.x = element_text(
                                             size = 7,
                                             angle = 45 + 90 + 90 + 90, hjust = .1,
                                             vjust = .9),
axis.text.y = element_text( size = 7,
                            angle = 0)) +
  scale_x_discrete(labels = fancy_axis) +
  scale_y_discrete(labels = fancy_axis)
ggsave(filename = "ACS AMI/pics/newp.jpeg", width = unit(5,"cm"), height = unit(5,"cm"))
