library(tidyverse)
rm(list = ls())
setwd('~/Desktop/ZhengPan/workspace/NUIST/John/胡文广/A4BX6/')
source('code/dependencies.R')


# txtfile <- read_lines('code/bandgap.txt')

getFeatureName <- function(txtfile){
  res7 <- lapply(tail(txtfile, 3), function(q){
    res1 <- str_replace_all(q, '\\"', '')
    res2 <- str_replace_all(res1, '\\[1\\]', '')
    res3 <- str_replace_all(res2, '\\[8\\]', '')
    res4 <- str_replace_all(res3, '\\[15\\]', '')
    res5 <- str_split(string = res4, pattern = ' ')[[1]]
    res6 <- res5[res5 != '']
  })
  res8 <- do.call(c, res7)
  return(res8)
}
getFeatureName(txtfile = txtfile)

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

changeName <- function(nameToChange){
  # nameToChange = formulaWithoutName[1]
  for (t in 20:1) {
    # t = 1
    pat <- str_c('feature_', t)
    nameToChange <- str_replace_all(string = nameToChange, 
                    pattern = pat, 
                    replacement = feaNames[t])  
  }
  return(nameToChange)
}


feaNames <- getFeatureName(txtfile = txtfile)

dotext <- function(df){
  # df= bandgapd
  # df= homod
  df <- as.data.frame(scale(df))
  formulaAndCorrelation <- lapply(1:length(txtfile), function(q){
    # q = 979
    if(str_detect(txtfile[q], 'Best Expression:')){
      formulaWithoutName <- strsplit(txtfile[q], split = 'Best Expression:  ')[[1]][2]
      formulaWithoutName <- strsplit(formulaWithoutName, ' ')[[1]]
      formulaWithName <- sapply(1:length(formulaWithoutName), function(s){
        changeName(formulaWithoutName[s])
      })
      res1 <- str_c(formulaWithName, collapse = '')
      
      attach(df)
      res2 <- evaluate::parse_all(res1)
      res3 <- eval(res2$expr[[1]])
      detach(df)
      res4 <- cor(df[, ncol(df)], res3)
      return(list(formula = res1, correlation = res4))
    }
  })
}

# aa <- dotext()

bandgap_formula <- sapply(1:length(aa), function(q){
  if(!is.null(aa[[q]])){
    return(aa[[q]]$formula)
  }
})
bandgap_formula <- do.call(c, bandgap_formula)
bandgap_cor <- sapply(1:length(aa), function(q){
  if(!is.null(aa[[q]])){
    return(aa[[q]]$correlation)
  }
})
bandgap_cor <- do.call(c, bandgap_cor)
bandgap <- data.frame(formula = bandgap_formula, correlation = bandgap_cor)
bandgap <- unique(bandgap)


final_step <- function(df, file_name){
  aa <- dotext(df = df)
  bandgap_formula <- sapply(1:length(aa), function(q){
    if(!is.null(aa[[q]])){
      return(aa[[q]]$formula)
    }
  })
  bandgap_formula <- do.call(c, bandgap_formula)
  bandgap_cor <- sapply(1:length(aa), function(q){
    if(!is.null(aa[[q]])){
      return(aa[[q]]$correlation)
    }
  })
  bandgap_cor <- do.call(c, bandgap_cor)
  bandgap <- data.frame(formula = bandgap_formula, correlation = bandgap_cor)
  bandgap <- unique(bandgap)
  # file_name = 'bandgap'
  write_csv(x = bandgap, file = str_c('results/', file_name, '.csv', collapse = ''))
  return(bandgap)
}

txtfile <- read_lines('code/homo.txt')
feaNames <- getFeatureName(txtfile = txtfile)
final_step(df = homod, file_name = 'homo')
txtfile <- read_lines('code/lumo.txt')
feaNames <- getFeatureName(txtfile = txtfile)
final_step(df = lumod, file_name = 'lumo')




