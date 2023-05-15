# preparation
setwd(dir = '~/Desktop/ZhengPan/workspace/NUIST/John/胡文广/A4BX6/')
library(tidyverse)
st1 <- read_csv(file = 'rawdata/A4BX6train_and_test_fill1.csv')
st2 <- st1[!is.na(st1$homo), ]
st3 <- st2[, c('A1_Symbol', 'A2_Symbol', 'A3_Symbol', 'A4_Symbol', 
               'B_Symbol', 
               'X1_Symbol', 'X2_Symbol', 'X3_Symbol', 'X4_Symbol', 'X5_Symbol', 'X6_Symbol')]
res4 <- sapply(1:nrow(st3), function(q){
  str_c(st3[q, ], collapse = '')
})
d1 <- lapply(1:10, function(q){
  return(dir(str_c('rawdata/update-0D/', q), full.names = T)  )
})
d1 <- do.call(c ,d1)
ds <- c(dir('rawdata/finish-update-0D', full.names = T), d1)

for(pat in ds[1:length(ds)]){
  path <- pat
  name_split <- getname(tail(strsplit(path, '/')[[1]], 1))
  if(str_c(name_split, sep = '', collapse = '') %in% res4){
    print(which(pat == ds)/length(ds)*100)
    path_of_folder <- str_c('rawdata/folders-2740', str_c(name_split, sep = '', collapse = ''), sep = '/', collapse = '')
    dir.create(path = path_of_folder)
     file.copy(from = str_c(pat, '/t.outmol'), to = path_of_folder)
  }
}
# write_csv(st1, file = 'rawdata/A4BX6train_and_test_fill.csv')
write_csv(st1, file = 'rawdata/A4BX6train_and_test_fill1.csv')
write_csv(st1, file = '~/Desktop/temp1.csv')

abbs <- sapply((colnames(st1)), function(q) {
  # q = colnames(st1)[1]
  sp1 <- strsplit(q, split = '')[[1]]
  res1 <- str_c(sp1[sp1 %in% c(LETTERS, as.character(1:10))], collapse = '')
  if(res1 == ''){
    res1 <- q
  }
  return(res1)
})

write_csv(data.frame(Abbreviation = abbs, Fullname = colnames(st1)), file = 'rawdata/Abbs.csv')








