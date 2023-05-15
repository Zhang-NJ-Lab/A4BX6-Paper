# preparation
st1 <- read_csv(file = 'rawdata/A4BX6train_and_test.csv')
st1$solvent = NA
st1$homo = NA
st1$lumo = NA
d1 <- lapply(1:10, function(q){
  return(dir(str_c('rawdata/update-0D/', q), full.names = T)  )
})
d1 <- do.call(c ,d1)
ds <- c(dir('rawdata/finish-update-0D', full.names = T), d1)

for(pat in ds[1:length(ds)]){
  # pat = ds[3869]
  # pat = ds[1]
  print(which(pat == ds))
  path <- pat
  name_split <- getname(tail(strsplit(path, '/')[[1]], 1))
  total_energy <- getinfo(path)
  if(is.na(total_energy)){
    total_energy <- list(NA, NA)
  }else{
    print(total_energy)
  }
  
  st1[st1$A1_Symbol == name_split["A1"] & st1$A2_Symbol == name_split["A2"] &
        st1$A3_Symbol == name_split["A3"] & st1$A4_Symbol == name_split["A4"] &
        st1$B_Symbol == name_split["B"] & st1$X1_Symbol == name_split["X1"] &
        st1$X2_Symbol == name_split["X2"] & st1$X3_Symbol == name_split["X3"] &
        st1$X4_Symbol == name_split["X4"] & st1$X5_Symbol == name_split["X5"] &
        st1$X6_Symbol == name_split["X6"], "solvent"] <- get_solvent_information(pat)
  st1[st1$A1_Symbol == name_split["A1"] & st1$A2_Symbol == name_split["A2"] &
       st1$A3_Symbol == name_split["A3"] & st1$A4_Symbol == name_split["A4"] &
       st1$B_Symbol == name_split["B"] & st1$X1_Symbol == name_split["X1"] &
       st1$X2_Symbol == name_split["X2"] & st1$X3_Symbol == name_split["X3"] &
       st1$X4_Symbol == name_split["X4"] & st1$X5_Symbol == name_split["X5"] &
       st1$X6_Symbol == name_split["X6"], "homo"] <- total_energy[[1]][1]
  st1[st1$A1_Symbol == name_split["A1"] & st1$A2_Symbol == name_split["A2"] &
       st1$A3_Symbol == name_split["A3"] & st1$A4_Symbol == name_split["A4"] &
       st1$B_Symbol == name_split["B"] & st1$X1_Symbol == name_split["X1"] &
       st1$X2_Symbol == name_split["X2"] & st1$X3_Symbol == name_split["X3"] &
       st1$X4_Symbol == name_split["X4"] & st1$X5_Symbol == name_split["X5"] &
       st1$X6_Symbol == name_split["X6"], "lumo"] <- total_energy[[2]][1]
  # st1 <- matchin(st1, pat)
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








