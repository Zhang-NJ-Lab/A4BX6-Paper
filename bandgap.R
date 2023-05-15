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
features <- st[, -ncol(st):(-ncol(st) + 1)] 
fea_new <- fea_slc_bycor(features, corr = 0.95)
homod <- cbind(fea_new, homo = st$homo)
lumod <- cbind(fea_new, lumo = st$lumo) 
bandgapd <- cbind(fea_new, bandgap = st$lumo - st$homo)
SR(df = bandgapd)