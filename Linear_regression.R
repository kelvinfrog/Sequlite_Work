library(wesanderson)
library(ggplot2)
library(plyr)
library(tidyr)
library(dplyr)
library(gtools)
library(broom)

#experiment number 
experiment = c("31k_Inc")
Channel = "Q3"
f_list = c("R3") #f_list = c("G1", "G2", "R3", "R4") <-use this if all four filters are needed
Cycle_number = 3


#file path for the tsv file and create path for the tsv file
file_path1 = 'C:/Users/kelvin/Documents/Bioinformatics_analysis/'
file_path2 = '_pixel_Q321.tsv'

experiment_list = sapply(experiment, toString)
file_list = sapply(experiment, function(x) paste0(file_path1, 'CRT', x, '/', 'CRT', x, file_path2), USE.NAMES=FALSE)
file_name_list = sapply(experiment, function(x) paste0('CRT', x), USE.NAMES=FALSE)


# Raw data fit
datalist = list()  

for (i in c(1:length(file_list)))
{
  #Read and reorganize data table
  db = read.table(file = file_list[i], 
                  sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  #Include this part if normalized data is needed
  db = db%>% filter(QC == "Q3") %>% group_by(Tile)
  db$Exp = file_name_list[i]
  datalist[[i]] <- db
}

db1 = do.call(rbind, datalist)

db1 = db1 %>% gather(Filter, FL, -Cycle,  -Tile, -QC, -Exp)
db1$Cycle <- gsub('Inc', '', db1$Cycle)
db1$Cycle <- as.integer(as.character(db1$Cycle))
db1$Tile <- gsub('mm', '', db1$Tile)
db1$Tile <- as.numeric(as.character(db1$Tile))


db2 = db1 %>% filter(QC == Channel) %>% group_by(Cycle, Exp, Filter) %>% summarise(mean_FL = mean(FL))


datalist2 = list()

# Fit all data
for (i in c(1 :length(f_list))) {
  db3 = db2 %>% filter(Filter == f_list[i], Cycle > Cycle_number)
  db4 = db3 %>% group_by(Exp)
  fitted_models =  db4 %>% do(model = lm(log(mean_FL) ~ Cycle, data = .)) 
  co_db = fitted_models %>% tidy(model) %>% filter(term == "Cycle") %>% select(Exp, estimate) %>% rename(coeff = estimate)
  r_db = fitted_models %>% glance(model) %>% select(Exp, r.squared)
  mdb =  merge(x = co_db, y = r_db, by = "Exp", all = TRUE)
  mdb$filter = f_list[i]  
  mdb$perdrop = (exp(mdb$coeff)-1)*100
  datalist2[[i]] <- mdb
}

#Combind the exponential fit data 
eb4 = do.call(rbind, datalist2)
eb4$model = "exponential"

datalist2 = list()
# Fit all data
for (i in c(1 :length(f_list))) {
  db3 = db2 %>% filter(Filter == f_list[i], Cycle > Cycle_number)
  db4 = db3 %>% group_by(Exp)
  fitted_models =  db4 %>% do(model = lm(mean_FL ~ Cycle, data = .)) 
  co_db = fitted_models %>% tidy(model) %>% filter(term == "Cycle") %>% select(Exp, estimate) %>% rename(coeff = estimate)
  r_db = fitted_models %>% glance(model) %>% select(Exp, r.squared)
  mdb =  merge(x = co_db, y = r_db, by = "Exp", all = TRUE)
  mdb$filter = f_list[i]
  mdb$perdrop = (exp(mdb$coeff)-1)*100
  datalist2[[i]] <- mdb
}

#Combind the exponential fit data 
lb4 = do.call(rbind, datalist2)
lb4$model = "linear"

datalist2 = list()
# Fit all data
for (i in c(1 :length(f_list))) {
  db3 = db2 %>% filter(Filter == f_list[i], Cycle > Cycle_number)
  db4 = db3 %>% group_by(Exp)
  fitted_models =  db4 %>% do(model = lm(mean_FL ~ log(Cycle), data = .)) 
  co_db = fitted_models %>% tidy(model) %>% filter(term == "log(Cycle)") %>% select(Exp, estimate) %>% rename(coeff = estimate)
  r_db = fitted_models %>% glance(model) %>% select(Exp, r.squared)
  mdb =  merge(x = co_db, y = r_db, by = "Exp", all = TRUE)
  mdb$filter = f_list[i]
  mdb$perdrop = (exp(mdb$coeff)-1)*100
  datalist2[[i]] <- mdb
}

#Combind the exponential fit data 
lob4 = do.call(rbind, datalist2)
lob4$model = "log"

total = rbind(eb4, lb4, lob4)
total$base.start = Cycle_number +1
View(total)








#Normalized fit
datalist = list()  

for (i in c(1:length(file_list)))
{
  #Read and reorganize data table
  db = read.table(file = file_list[i], 
                  sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  #Include this part if normalized data is needed
  db = db%>%
    filter(QC == "Q3") %>% group_by(Tile) %>%
    #mutate(G1n = G1/G1[Cycle =='Inc1']) %>%
    mutate(G2n = G2/G2[Cycle =='Inc1']) %>%
    #mutate(R3n = R3/R3[Cycle =='Inc1']) %>%
    #mutate(R4n = R4/R4[Cycle =='Inc1']) %>%
    select(Cycle, Tile,  QC, G2n) %>%
    #select(Cycle, Tile,  QC, G1n, G2n, R3n, R4n) %>%
    ungroup(Tile) %>%
    #rename(G1 = G1n, G2=G2n, R3=R3n, R4=R4n)
    rename(G2 = G2n)
  
  db$Exp = file_name_list[i]
  datalist[[i]] <- db
}

db1 = do.call(rbind, datalist)

db1 = db1 %>% gather(Filter, FL, -Cycle,  -Tile, -QC, -Exp)
db1$Cycle <- gsub('Inc', '', db1$Cycle)
db1$Cycle <- as.integer(as.character(db1$Cycle))
db1$Tile <- gsub('mm', '', db1$Tile)
db1$Tile <- as.numeric(as.character(db1$Tile))


db2 = db1 %>% filter(QC == Channel) %>% group_by(Cycle, Exp, Filter) %>% summarise(mean_FL = mean(FL))

datalist2 = list()

# Fit all data
for (i in c(1 :length(f_list))) {
  db3 = db2 %>% filter(Filter == f_list[i], Cycle > Cycle_number)
  db4 = db3 %>% group_by(Exp)
  fitted_models =  db4 %>% do(model = lm(log(mean_FL) ~ Cycle, data = .)) 
  co_db = fitted_models %>% tidy(model) %>% filter(term == "Cycle") %>% select(Exp, estimate) %>% rename(coeff = estimate)
  r_db = fitted_models %>% glance(model) %>% select(Exp, r.squared)
  mdb =  merge(x = co_db, y = r_db, by = "Exp", all = TRUE)
  mdb$filter = f_list[i]  
  datalist2[[i]] <- mdb
}

#Combind the exponential fit data 
eb4 = do.call(rbind, datalist2)
eb4$model = "exponential"

datalist2 = list()
# Fit all data
for (i in c(1 :length(f_list))) {
  db3 = db2 %>% filter(Filter == f_list[i], Cycle > Cycle_number)
  db4 = db3 %>% group_by(Exp)
  fitted_models =  db4 %>% do(model = lm(mean_FL ~ Cycle, data = .)) 
  co_db = fitted_models %>% tidy(model) %>% filter(term == "Cycle") %>% select(Exp, estimate) %>% rename(coeff = estimate)
  r_db = fitted_models %>% glance(model) %>% select(Exp, r.squared)
  mdb =  merge(x = co_db, y = r_db, by = "Exp", all = TRUE)
  mdb$filter = f_list[i]  
  datalist2[[i]] <- mdb
}

#Combind the exponential fit data 
lb4 = do.call(rbind, datalist2)
lb4$model = "linear"

datalist2 = list()
# Fit all data
for (i in c(1 :length(f_list))) {
  db3 = db2 %>% filter(Filter == f_list[i], Cycle > Cycle_number)
  db4 = db3 %>% group_by(Exp)
  fitted_models =  db4 %>% do(model = lm(mean_FL ~ log(Cycle), data = .)) 
  co_db = fitted_models %>% tidy(model) %>% filter(term == "log(Cycle)") %>% select(Exp, estimate) %>% rename(coeff = estimate)
  r_db = fitted_models %>% glance(model) %>% select(Exp, r.squared)
  mdb =  merge(x = co_db, y = r_db, by = "Exp", all = TRUE)
  mdb$filter = f_list[i]  
  datalist2[[i]] <- mdb
}

#Combind the exponential fit data 
lob4 = do.call(rbind, datalist2)
lob4$model = "log"

total = rbind(eb4, lb4, lob4)
total$base.start = Cycle_number +1
View(total)

