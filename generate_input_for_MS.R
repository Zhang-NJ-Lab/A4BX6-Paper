library(tidyverse)
library(reshape2)
library(readr)
rm(list = ls())
setwd("~/Desktop/ZhengPan/workspace/John/胡文广/MS文件生成/")
#### folders construction  =================================================================
sol <- as.data.frame(read_csv(file = "Rawdata/Solvant.csv",col_names = FALSE))
colnames(sol) <- c('value', 'name')
sol$value <- format(sol$value, nsmall = 4)
car1 <- read_lines(file = "Rawdata/p416.car")
input1 <- read_lines(file = "Rawdata/p416.input")
sol_line <- sapply(1:nrow(sol), function(q){
  temp <- strsplit(input1[31], ' ')[[1]]
  temp[15] <- sol$value[q]
  temp[20] <- sol$name[q]
  str_c(temp, collapse = ' ')
})

# input1[31] <- sol_line[6]

# write_lines(car1, file = 't.car')
# write_lines(input1, file = 't.input')

change_car <- function(f_element, f_car_to_c){
  # f_car_to_c = car1
  # f_element = rep("AA", 11)
  # q = 100
  # f_element = strsplit(sel_A4BX6[q], ',')[[1]]
  # f_car_to_c = car1
  # a <- change_car(, car1)
  f_element[c(1,6,length(f_element), length(f_element) - 1, 2, 3, 4, 5, 7, 8, 9)] <- f_element
  line_index <- which(str_length(f_car_to_c) > 30)
  if(length(f_element) != length(line_index)){
    print("incorrect length of input")
    return(NULL)
  }
  for (f_q in line_index) {
    # f_q = 12
    # f_q = 5
    q_str <- f_car_to_c[f_q]
    q_split1 <- strsplit(q_str, ' ')[[1]]
    q_split1 <- q_split1[q_split1 != '']
    q_split1[length(q_split1) - 1] <- f_element[which(f_q == line_index)]
    q_split1[1] <- paste0(f_element[which(f_q == line_index)], which(f_q == line_index))
    if(str_length(q_split1[length(q_split1) - 1]) == 1){
      q_split1[length(q_split1) - 1] <- str_c(q_split1[length(q_split1) - 1], ' ', collapse = '')
    }
    q_split1[1] <- str_c(q_split1[1], str_c(rep(' ', 4 - str_length(q_split1[1])), collapse = ''), collapse = '')
    # f_car_to_c[f_q] <- str_c(q_split1, sep = '', collapse = ' ')
    
    if(f_q %in% c(8, 12, 13, 14, 15)){
      f_car_to_c[f_q] <- paste0(q_split1[1], '    ', q_split1[2], '   ', q_split1[3], '    ',
                                q_split1[4], ' ', q_split1[5], ' ', q_split1[6], '      ',q_split1[7], '      ',
                                q_split1[8], '  ',q_split1[9])
    }else{
      f_car_to_c[f_q] <- paste0(q_split1[1], '    ', q_split1[2], '   ', q_split1[3], '   ',
                                q_split1[4], ' ', q_split1[5], ' ', q_split1[6], '      ',q_split1[7], '      ',
                                q_split1[8], '  ',q_split1[9])
    }
  }
  return(f_car_to_c)
}
# a <- change_car(rep("AA", 11), car1)
# write_lines(a, file = 'Results/a.car')

#===== construct A4BX6 ####
# A: Li Na K Rb Cs
# B: Ge Sn Pb
# X: F Cl Br I

od <- function(f_s){
  # f_s = "3,2,1,1"
  # str_c(sort((as.numeric(strsplit(f_s, ',')[[1]]))), collapse = ',')
}
conver_ele <- function(f_s,f_pos){
  # f_s = "1,1,2,3"
  # f_pos = "A"
  if(f_pos == "A"){
    str_c(A[as.numeric(strsplit(f_s, ',')[[1]])], collapse = ',')
  }else if(f_pos == "X"){
    str_c(X[as.numeric(strsplit(f_s, ',')[[1]])], collapse = ',')
  }else{
    NULL
  }
}

A <- c('Li', 'Na', 'K', 'Rb', 'Cs')
B <- c('Ge', 'Sn', 'Pb')
X <- c('F', 'Cl', 'Br', 'I')

A_arr <- array(rep(0, 5^4), c(5,5,5,5))
for (s1 in 1:5) {
  for (s2 in 1:5) {
    for (s3 in 1:5) {
      for (s4 in 1:5) {
            A_arr[s1, s2, s3, s4] <- od(paste(s1, s2, s3, s4, sep = ','))
      } 
    } 
  }
}
A_vec <- as.vector(A_arr)
A_vec <- A_vec[!duplicated(A_vec)]
A_vec <- sapply(1:length(A_vec), function(q){
  conver_ele(A_vec[q], 'A')
})

X_arr <- array(rep(0, 4^6), c(4,4,4,4,4,4))
for (s1 in 1:4) {
  for (s2 in 1:4) {
    for (s3 in 1:4) {
      for (s4 in 1:4) {
        for (s5 in 1:4) {
          for (s6 in 1:4) {
            X_arr[s1, s2, s3, s4, s5, s6] <- od(paste(s1, s2, s3, s4, s5, s6, sep = ','))
          } 
        } 
      } 
    } 
  }
}
X_vec <- as.vector(X_arr)
X_vec <- X_vec[!duplicated(X_vec)]
X_vec <- sapply(1:length(X_vec), function(q){
  # q = 1
  conver_ele(X_vec[q], 'X')
})

A4BX6 <- array(0, c(length(A_vec), length(B), length(X_vec)))
for (a in 1:length(A_vec)) {
  for (b in 1:length(B)) {
    for (x in 1:length(X_vec)) {
      A4BX6[a, b, x] <-  str_c(A_vec[a], B[b], X_vec[x], sep = ',')
    }
  }
}
A4BX6 <- as.vector(A4BX6)
length(A4BX6)
#=====

ind <- round(runif(10000, min = 0, max = 1) * (length(A4BX6) - 1), 0) + 1
ind <- ind[!duplicated(ind)][1:5000]
sel_A4BX6 <- A4BX6[ind]

for(q in 1:length(sel_A4BX6)){
  # q = 100
  print(paste0('current job: ', q, '  ', q/length(sel_A4BX6)*100, ' %'))
  a <- change_car(strsplit(sel_A4BX6[q], ',')[[1]], car1)
  l_name <- str_c(strsplit(sel_A4BX6[q], ',')[[1]], collapse = '')
  setwd("~/Desktop/ZhengPan/workspace/John/MS文件生成/Results/")
  folder_name <- paste0(l_name, '/')
  system(paste0('mkdir ', folder_name))
  setwd(paste0("~/Desktop/ZhengPan/workspace/John/MS文件生成/Results/", folder_name))
  write_lines(a, file = paste0('t.car'))
  input1[31] <- sol_line[round(runif(1,min = 1, max = 15), 0)]
  write_lines(input1, file = 't.input')
}
setwd("~/Desktop/ZhengPan/workspace/John/MS文件生成/")
system('zip -r dmol3_files.zip Results/')

#===== train and test =======
raw_ele <- read_csv(file = 'Rawdata/Periodic Table of Elements.csv')
folder_names <- dir(path = 'Results/folders/')
extract_name <- function(f_names){
  # f_names = folder_names[1]
  sp <- strsplit(f_names, split = '')[[1]]
  res1 <- sapply(which(sp %in% LETTERS), function(q){
    # q = 1
    if(q == length(sp)){
      return(sp[q])
    }else if(sp[q + 1] %in% letters){
      return(str_c(sp[q], sp[q + 1], collapse = ''))
    }else{
      return(sp[q])
    }
  })
}
A4BX6_string <- lapply(1:length(folder_names), function(q){
  extract_name(folder_names[q])
})
A4BX6_string <- do.call(rbind, A4BX6_string)
colnames(A4BX6_string) <- c('A1', 'A2', 'A3', 'A4', 'B', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6')
A4BX6_string <- as.data.frame(A4BX6_string)
train_and_test_l <- lapply(1:nrow(A4BX6_string), function(q){
  # q = 1
  print(paste0('current job: ', q, '   ', round(q/nrow(A4BX6_string), 3)*100, ' % '))
  t1 <- A4BX6_string[q, ]
  res1 <- lapply(1:ncol(t1), function(s){
    # s = 2
    res <- filter(raw_ele, Symbol == t1[, s]) %>%
      select(-(AtomicNumber:Element))
    if(s <= 4){
      colnames(res) <- str_c(paste0('A', s ,'_'), colnames(res), sep = '')
    }else if(s == 5){
      colnames(res) <- str_c(paste0('B','_'), colnames(res), sep = '')
    }else{
      colnames(res) <- str_c(paste0('X', s - 5 ,'_'), colnames(res), sep = '')
    }
    return(as.data.frame(res))
  })
  res1 <- do.call(cbind, res1)
})
train_and_test <- do.call(rbind, train_and_test_l)
filter(raw_ele, Symbol == 'He')
write_csv(train_and_test, file = 'results/train_and_test.csv')

train_and_test <- read_csv(file = 'results/train_and_test.csv')
#========================================================================

#### 

#### test ####


