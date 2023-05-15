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
# ggsave(filename = "results/feature_engineering/feature_ranking_homo.jpeg", plot = cor_for_fea_sel_df_pic_homo, width = unit(7,"in"), height = unit(10,"in"))
# ggsave(filename = "results/feature_engineering/feature_ranking_lumo.jpeg", plot = cor_for_fea_sel_df_pic_lumo, width = unit(7,"in"), height = unit(10,"in"))

