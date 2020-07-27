library(wesanderson)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(stringr)
library(gtools)
library(ggpmisc)
library(broom)
library(grid)
library(gridExtra)
library(qdapRegex)
library(gtable)
library(purrr)
library(data.table)


# Input -------------------------------------------------------------------

#experiment number, only the number with NO CRT
experiment  = c("106x_Inc", "106x_CL", "107x_Inc", "107x_CL", "108x_Inc", "36y_CL", "36y_Inc")

file_path1 = 'C:/Users/kelvin/Documents/Bioinformatics_analysis/'
file_path2 = 'C:/Users/kelvin/Documents/Figures/'

exp = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1), pattern = e, full.names = TRUE, recursive = F)
  exp = c(exp, x)
}

experiment = gsub(file_path1, "", exp)


x.font = 15
y.font = 15
theme.font = 13

numrow = function(df){
  x = length(unique(df$Tile))
  if (x <= 5){
    row = 1
  } else if (x > 5 & x <=10) {
    row = 2
  } else if (x > 10 & x <=15) {
    row = 3 
  } else if (x > 15 & x <= 20) {
    row = 4 
  } else if (x > 20) {
    row = 5 
  }
  return(row)
}


# SNR analysis ------------------------------------------------------------

#file path for the tsv file and create path for the tsv file
SNR_path = '_fbqc_SNRMean.tsv'


SNR_list = sapply(experiment, function(x) paste0(file_path1,  x, '/fbqc/',  x, SNR_path), USE.NAMES=FALSE)
file_name_list = sapply(experiment, function(x) paste0( x), USE.NAMES=FALSE)

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (e in experiment){
  f_list = paste0(file_path2, 'SNR_plot/',  e, '_SNR_plot/')
  dir_list = c(dir_list, f_list)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)


SNR_table = list() 

for (i in c(1:length(SNR_list))){
  db1 = read.table(file = SNR_list[i], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db1$Exp = file_name_list[i]
  SNR_table[[i]] <- db1
}

sdb = do.call(rbind, SNR_table)

db1 =  sdb
db1$Cycle <- gsub('Inc', '', db1$Cycle)
db1$Cycle <- as.integer(as.character(db1$Cycle))
db1 = db1 %>% select(-Tile.1)

dbm = db1  %>% group_by(Exp, Base, Cycle)  %>% summarise(QC.Value = mean(QC.Value))
dbm$Tile = "Mean"
dbme = db1  %>% group_by(Exp, Base, Cycle)  %>% summarise(QC.Value = median(QC.Value))
dbme$Tile = "Median"
db2 =  bind_rows(db1, dbm, dbme)
db2$Tile <- factor(db2$Tile, levels=mixedsort(unique(db2$Tile)))

db2$Base[db2$Base == "A"] <- "G1"
db2$Base[db2$Base == "T"] <- "G2"
db2$Base[db2$Base == "G"] <- "R3"
db2$Base[db2$Base == "C"] <- "R4"
db2$Tile <- factor(db2$Tile, levels=mixedsort(unique(db2$Tile)))


for (i in c(1 :length(file_name_list))){
  
  #Wrap plot
  db3 = db2 %>% filter(Exp == file_name_list[i])
  num.row = numrow(db3)
  
  P1 = ggplot(db2 %>% filter(Exp == file_name_list[i]), aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 0.7, alpha = 0.7) + 
    geom_point(size=0.5, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size= theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Signal to Noise Ratio"), 
         x = "Cycle", y = "SNR") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~. , nrow = num.row)
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'SNR', '_wrap_plot', '.png'),
         width = 40, height = 20, units = "cm", plot = P1)
  
  
  
  #Mean plot
  P2 = ggplot(db2 %>% filter(Exp == file_name_list[i], Tile == "Mean"), aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=2, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text( size = 20)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = 18)) +
    theme(strip.text = element_text(size=20)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Signal to Noise Ratio"), 
         x = "Cycle", y = "SNR") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~. , ncol = 1)
  
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'Mean_SNR', '.png'),
         width = 40, height = 20, units = "cm", plot = P2)
  
  
}

# RFL vs Noise analysis ---------------------------------------------------

#file path for the tsv file and create path for the tsv file
relative_path = '_fbqc_RelativeIntensityMean.tsv'
noise_path = '_fbqc_NoiseMean.tsv'


relative_list = c()
noise_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = '_fbqc_RelativeIntensityMean.tsv$', full.names = TRUE, recursive = TRUE)
  y = list.files(path =  paste0(file_path1,  e),
                 pattern = '_fbqc_NoiseMean.tsv$', full.names = TRUE, recursive = TRUE)
  relative_list = c(relative_list, x)
  noise_list = c(noise_list, y)
}

file_name = unlist(rm_between(relative_list, 'Bioinformatics_analysis/', '/', extract=TRUE))

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (i in c(1: length(file_name))){
  x = paste0('C:/Users/kelvin/Documents/Figures/SNR_plot/', file_name[i], '_SNR_plot/')
  dir_list = c(dir_list, x)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)


for (f in c(1 : length(relative_list))){
  dfr = read.delim(file = relative_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  dfn = read.delim(file = noise_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  dfr$QC = 'RFL'
  dfn$QC = 'Noise'
  
  df = rbind(dfr, dfn)
  df$Cycle = gsub("Inc", "", df$Cycle)
  df = df %>%  group_by(Cycle, Base, QC) %>% summarise(Value = mean(QC.Value))
  df$Cycle = as.numeric(df$Cycle)
  df$Base[df$Base == "A"] <- "G1"
  df$Base[df$Base == "T"] <- "G2"
  df$Base[df$Base == "G"] <- "R3"
  df$Base[df$Base == "C"] <- "R4"
  df$group = paste0(df$Base, "_", df$QC)
  
  datalist1 = list()
  channel = c("G1", "G2", "R3", "R4")
  for (i in c(1 :length(channel))) {
    r.vector = c()
    db1 = df %>% filter(Base == channel[i], Value != 0) %>% group_by(QC)
    fitted_models =  db1 %>% do(model = lm(log(Value) ~ Cycle, data = .)) 
    co_db = fitted_models %>% ungroup() %>% transmute(QC, Coef = map(model, tidy)) %>% unnest(Coef) %>% 
      filter(term == "Cycle") %>% select(QC, estimate)
    for (c in c(1: nrow(fitted_models))) r.vector[c] <- summary(fitted_models$model[[c]])$r.squared
    co_db$r.squared = r.vector
    co_db$filter = channel[i]  
    datalist1[[channel[i]]] <- co_db
  } 
  
  
  #Combind the exponential fit data 
  db2 = bind_rows( datalist1)
  db2$percent_change = (exp(db2$estimate)-1)
  db2$group = paste0(db2$filter, "_", db2$QC)
  db3 = db2 %>% select(percent_change, group, r.squared)
  db3$percent_change = round(db3$percent_change *100, 2)
  db3$r.squared = round(db3$r.squared, 3)
  db4 = as.data.frame(t(as.matrix(db3)))
  colnames(db4) = db3$group
  db4 = db4[c(1,3),]
  db4 = setDT(db4, keep.rownames = TRUE)[]
  colnames(db4) = c('', db3$group)
  
  tt3 <- ttheme_minimal(
    core=list(bg_params = list(blues9[1:6], col=NA),
              fg_params=list(fontface=2, fontsize = 15)),
    colhead=list(fg_params=list(col="navyblue", fontface=4, fontsize = 15)))
  
  tbl =  tableGrob(db4, rows =  NULL, theme = tt3)
  
  
  #This is for plotting
  df = df[complete.cases(df), ]
  di = mean(df$Value[df$QC == "Noise"])/mean(df$Value[df$QC == "RFL"])
  
  
  
  df1 = df %>% spread(QC, Value)
  df1$Noise1 = df1$Noise/di
  df2 = df1 %>% select(-Noise) %>% gather("QC", "Value", - Base, -Cycle, -group) %>% drop_na()
  df2$QC[df2$QC == "Noise1"] <- "Noise"
  
  x.font=15
  y.font=15
  theme.font = 18
  
  P1 = ggplot(df2, aes(x= Cycle, y = Value, col = QC)) +
    theme_bw() +
    geom_point(size=1.5, alpha = 0.7) + 
    geom_smooth(method="glm",
                method.args=list(family=gaussian(link="log")), se = FALSE) +
    scale_y_continuous(sec.axis = sec_axis(~.* di, name="Noise")) + 
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size = theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name[f], " RFL vs. Noise"), x = "Cycle", y = "RFL") +
    #theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Base ~., ncol = 2, scales = "free_y")
  
  g <- ggplot_gtable(ggplot_build(P1))
  stript <- which(grepl('strip-t', g$layout$name))
  
  n = length(unique(db5$base))
  fills <- c( rep("#4363d8", 1), rep("#f032e6", 1), rep("#3cb44b", 1), rep("#f58231", 1)) 
  
  k <- 1
  
  for (i in stript) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  
  #P2 = grid.draw(g)
  
  P3a = grid.arrange(tbl, P1, 
                     nrow=2, 
                     heights=c(2, 8))
  
  ggsave(filename = paste0(dir_list[f], file_name[f], '_', '_RFL_vs_noise','.png'), 
         width = 40, height = 20, units = "cm", plot = P3a)
}


# FOG Analysis ------------------------------------------------------------

#file path for the tsv file and create path for the tsv file
FOG_path = '_fbqc_Fractionofgoodfeatures.tsv'


FOG_list = sapply(experiment, function(x) paste0(file_path1,  x, '/fbqc/',  x, FOG_path), USE.NAMES=FALSE)
file_name_list = sapply(experiment, function(x) paste0( x), USE.NAMES=FALSE)

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (e in experiment){
  f_list = paste0(file_path2, 'FOG_plot/',  e, '_FOG_plot/')
  dir_list = c(dir_list, f_list)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)


FOG_table = list() 

for (i in c(1:length(FOG_list))){
  db1 = read.table(file = FOG_list[i], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db1$Exp = file_name_list[i]
  FOG_table[[i]] <- db1
}

sdb = do.call(rbind, FOG_table)

db1 =  sdb
db1$Cycle <- gsub('Inc', '', db1$Cycle)
db1$Cycle <- as.integer(as.character(db1$Cycle))
db1 = db1 %>% select(-Tile.1)

dbm = db1  %>% group_by(Exp, Base, Cycle)  %>% summarise(QC.Value = mean(QC.Value))
dbm$Tile = "Mean"
dbme = db1  %>% group_by(Exp, Base, Cycle)  %>% summarise(QC.Value = median(QC.Value))
dbme$Tile = "Median"
db2 =  bind_rows(db1, dbm, dbme)
db2$Tile <- factor(db2$Tile, levels=mixedsort(unique(db2$Tile)))

db2$Base[db2$Base == "A"] <- "G1"
db2$Base[db2$Base == "T"] <- "G2"
db2$Base[db2$Base == "G"] <- "R3"
db2$Base[db2$Base == "C"] <- "R4"
#db2 <- rename(db2, FOG = QC.Value)
db2$Tile <- factor(db2$Tile, levels=mixedsort(unique(db2$Tile)))


for (i in c(1 :length(file_name_list))){
  
  #Wrap plot
  db3 = db2 %>% filter(Exp == file_name_list[i])
  num.row = numrow(db3)
  
  P1 = ggplot(db3 %>% filter(Exp == file_name_list[i]), aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 1, alpha = 0.7) + 
    geom_point(size=1, alpha = 0.7) + 
    scale_y_continuous(labels=scales::percent) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size= x.font)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size= y.font)) +
    theme(strip.text = element_text(size= theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Fraction of Good Features"), 
         x = "Cycle", y = "Fraction of Good Features") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~. , nrow = num.row)
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'FOG', '_wrap_plot', '.png'),
         width = 40, height = 20, units = "cm", plot = P1)
  
  
  #Mean plot
  P2 = ggplot(db2 %>% filter(Exp == file_name_list[i], Tile == "Mean"), aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=2, alpha = 0.7) + 
    scale_y_continuous(labels=scales::percent) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text( size = 20)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = 18)) +
    theme(strip.text = element_text(size=20)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Fraction of Good Features"), 
         x = "Cycle", y = "Fraction of Good Features") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~. , ncol = 1)
  
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'Mean_FOG', '.png'),
         width = 40, height = 20, units = "cm", plot = P2)
  
}


# Pixel Analysis ----------------------------------------------------------


#file path for the tsv file and create path for the tsv file
pixel_path = '_pixel_percentile_group.tsv'


#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (e in experiment){
  f_list = paste0(file_path2, 'Pixel_plot/',  e, '_Pixel_plot/')
  dir_list = c(dir_list, f_list)
}

#experiment_list = sapply(experiment, toString)
file_list = sapply(experiment, function(x) paste0(file_path1,  x, '/fbqc/',  x, pixel_path), USE.NAMES=FALSE)
file_name_list = sapply(experiment, function(x) paste0( x), USE.NAMES=FALSE)

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)

Channel = "Q3"
f_list = c("G1", "G2", "R3", "R4")

datalist = list()  

for (i in c(1:length(file_list))){
  db = read.delim(file = file_list[i], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db = na.omit(db)
  db$Exp = file_name_list[i]
  
  db = db %>% 
    select(Cycle, Color, Tile, Exp, X25th, X50th, X75th) %>% 
    gather(QC, FL, -Cycle, -Color, -Tile, -Exp)
  colnames(db) = c("Cycle", "Filter", "Tile",  "Exp",   "QC", "FL")
  db$Tile[db$Tile == "0mm"] <- "00mm"
  db$Tile[db$Tile == "5mm"] <- "05mm"
  db$QC[db$QC == "X25th"] <- "Q1"
  db$QC[db$QC == "X50th"] <- "Q2"
  db$QC[db$QC == "X75th"] <- "Q3"
  db = spread(db, Filter, FL)
  
  db = db %>% gather(Filter, FL, -Cycle, -QC, -Tile, -Exp)
  db$Cycle <- gsub('Inc', '', db$Cycle)
  db$Cycle <- as.integer(as.character(db$Cycle))
  datalist[[i]] <- db
}

db1 = do.call(rbind, datalist)

datalist2 = list()

experiment_list = unique(db1$Exp)

for (e in c(1 :length(experiment_list))) {
  for (i in c(1 :length(f_list))) {
    r.vector = c()
    db2 = db1 %>% filter(QC == Channel, Filter == f_list[i], Exp == experiment_list[e])
    db3 = db2 %>% group_by(Tile)
    fitted_models =  db3 %>% do(model = lm(log(FL) ~ Cycle, data = .)) 
    co_db = fitted_models %>% ungroup() %>% transmute(Tile, Coef = map(model, tidy)) %>% unnest(Coef) %>% 
      filter(term == "Cycle") %>% select(Tile, estimate)
    for (c in c(1: nrow(fitted_models))) r.vector[c] <- summary(fitted_models$model[[c]])$r.squared
    co_db$r.squared = r.vector
    co_db$filter = f_list[i]  
    datalist2[[paste0(f_list[i],"_", experiment_list[e])]] <- co_db
  }
}


#Combind the exponential fit data 
db4 = bind_rows( datalist2, .id = 'Exp')
db4$perchange = (exp(db4$estimate)-1)

#Create a column to specifify the experiment name
exp_col = gsub("G\\d_", "", db4$Exp)
exp_col = gsub("R\\d_", "", exp_col)
exp_col = gsub("\\.\\d+", "", exp_col)
db4$Exp = exp_col
db4$Tile = gsub("mm", "", db4$Tile)
#db4$Tile = as.numeric(db4$Tile)

#add mean and median data to the plot
db5 = db4 %>% group_by(filter, Exp) %>% summarize(perchange = mean(perchange), estimate = mean(estimate))
db5b = db4 %>% group_by(filter, Exp) %>% summarize(perchange = median(perchange), estimate = median(estimate))
db5$Tile = "mean"
db5b$Tile = "median"
db4$Tile = as.character(db4$Tile)
db6 =  bind_rows(db5, db5b, db4 %>% select(-r.squared))
db6$Tile <- factor(db6$Tile, levels=mixedsort(unique(db6$Tile)))


#Print Pixel figure Cycle vs, FL
for (i in c(1:length(file_name_list))){
  
  dbm = db1 %>% filter(QC == "Q3", Exp == file_name_list[i] ) %>% group_by(Filter, Cycle, QC, Exp)  %>% summarise(FL = mean(FL))
  dbm$Tile = "mean"
  dbme = db1 %>% filter(QC == "Q3", Exp == file_name_list[i]) %>% group_by(Filter, Cycle, QC, Exp)  %>% summarise(FL = median(FL))
  dbme$Tile = "median"
  db1a =  bind_rows(db1 %>% filter(QC == "Q3", Exp == file_name_list[i]), dbm, dbme)
  db1a$Tile <- factor(db1a$Tile, levels=mixedsort(unique(db1a$Tile)))
  num.row = numrow(db1a)
  
  P1 = ggplot(db1a %>% filter(QC == "Q3"), aes(x = Cycle, y = FL, group = Filter, colour=Filter)) + 
    theme_bw() +
    geom_line(size = 1, alpha = 0.7) + 
    geom_point(size=1, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size = theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0( experiment[i], " Pixel Intensity (75% percentile)"), 
         x = "Cycle", y = "Pixel Intensity") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~., nrow = num.row)
  
  
  ggsave(filename = paste0(dir_list[i], file_name_list[i], '_', "Pixel_Intensity_Six_Tile", '_wrap_plot', '.png'),
         width = 40, height = 20, units = "cm", plot = P1)
  
  P1a = ggplot(db1a %>% filter(QC == "Q3"), aes(x = Cycle, y = FL, group = Filter, colour=Filter)) + 
    theme_bw() +
    geom_line(size = 1, alpha = 0.7) + 
    geom_point(size=1, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size = theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0( experiment[i], " Pixel Intensity (75% percentile)"), 
         x = "Cycle", y = "Pixel Intensity") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~., nrow = num.row, scales = "free_y")
  
  
  ggsave(filename = paste0(dir_list[i], file_name_list[i], '_', "Pixel_Intensity_Six_Tile", '_wrap_plot_2', '.png'),
         width = 40, height = 20, units = "cm", plot = P1a)
  
  P2 = ggplot(db1a %>% filter(QC == "Q3", Tile == "mean"), aes(x = Cycle, y = FL, group = Filter, colour=Filter)) + 
    theme_bw() +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=2, alpha = 0.7) + 
    #geom_smooth(method = "lm", se=FALSE, formula = log(y) ~ x,  alpha = 0.3) +
    # stat_poly_eq(formula = log(y) ~ x, 
    #              aes(label = paste(..eq.label.., sep = "~~~")), 
    #              label.x = "left", label.y = "top", parse = TRUE, size = 6) +  
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text( size = 20)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = 20)) +
    theme(strip.text = element_text(size=28)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0( experiment[i], " Pixel Intensity (75% percentile)"), 
         x = "Cycle", y = "Pixel Intensity") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~., ncol = 3)
  
  
  ggsave(filename = paste0(dir_list[i], file_name_list[i], '_', 'Pixel_Intensity_Mean', '.png'),
         width = 40, height = 20, units = "cm", plot = P2)
  
  P3 = ggplot(db6 %>% filter(Exp == file_name_list[i]), aes(x = Tile, y = perchange, colour=filter, group = filter)) + 
    theme_bw() +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=4, alpha = 0.7) + 
    scale_y_continuous(labels=scales::percent) +
    labs(title = "Pixel intensity percentage change per cycle", x = "Tile", y = "Percentage change per cycle") +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6"), 
                       name = "Filter") +
    theme(plot.title = element_text(hjust = 0.5, size=25)) +
    theme(axis.title.x=element_text(size=25), axis.text.x = element_text(size = 18, angle = 45, hjust = 0.9)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=25), 
          axis.text.y = element_text(size = 22)) +
    theme(legend.title = element_text(size = 25), legend.text = element_text(size = 22)) +
    theme(strip.text = element_text(size=25)) +
    facet_grid(. ~ Exp)
  
  db7 = db6 %>% filter(Exp == file_name_list[i]) %>% select(filter, perchange, Tile)
  names(db7) = c("Filter", "Percent_change", "Tile")
  db7 = db7 %>% spread(Tile, Percent_change)
  db7[, -1] = db7[, -1]*100
  db7 <- db7 %>% 
    mutate_if(is.numeric, round, digits = 2)
  
  tt3 <- ttheme_minimal(
    core=list(bg_params = list(fill = c("#3cb44b", "#f58231", "#4363d8", "#f032e6"), col=NA),
              fg_params=list(fontface=3)),
    colhead=list(fg_params=list(col="navyblue", fontface=4)))
  
  tbl =  tableGrob(db7, rows =  NULL, theme = tt3)
  
  P3a = grid.arrange(P3, tbl, 
                     nrow=2, 
                     heights=c(8, 2))
  
  ggsave(filename = paste0(dir_list[i], file_name_list[i], '_', 'Pixel_Percent_Change_1', '.png'),
         width = 40, height = 20, units = "cm", plot = P3a)
  
  P4 = ggplot(db6 %>% filter(Exp == file_name_list[i], !Tile %in% c('mean', 'median') ), aes(x = Tile, y = perchange, colour=filter, group = filter)) + 
    theme_bw() +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=4, alpha = 0.7) + 
    scale_y_continuous(labels=scales::percent) +
    labs(title = "Pixel intensity percentage change per cycle", x = "Tile", y = "Percentage change per cycle") +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6"), 
                       name = "Filter") +
    theme(plot.title = element_text(hjust = 0.5, size=25)) +
    theme(axis.title.x=element_text(size=25), axis.text.x = element_text(size = 18, angle = 45, hjust = 0.9)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=25), 
          axis.text.y = element_text(size = 22)) +
    theme(legend.title = element_text(size = 25), legend.text = element_text(size = 22)) +
    theme(strip.text = element_text(size=25)) +
    facet_grid(. ~ Exp)
  
  db7 = db6 %>% filter(Exp == file_name_list[i]) %>% select(filter, perchange, Tile)
  names(db7) = c("Filter", "Percent_change", "Tile")
  db7 = db7 %>% spread(Tile, Percent_change)
  db7[, -1] = db7[, -1]*100
  db7 <- db7 %>% 
    mutate_if(is.numeric, round, digits = 2)
  
  tt3 <- ttheme_minimal(
    core=list(bg_params = list(fill = c("#3cb44b", "#f58231", "#4363d8", "#f032e6"), col=NA),
              fg_params=list(fontface=3)),
    colhead=list(fg_params=list(col="navyblue", fontface=4)))
  
  tbl =  tableGrob(db7, rows =  NULL, theme = tt3)
  
  P4a = grid.arrange(P4, tbl, 
                     nrow=2, 
                     heights=c(8, 2))
  
  ggsave(filename = paste0(dir_list[i], file_name_list[i], '_', 'Pixel_Percent_Change_2', '.png'),
         width = 40, height = 20, units = "cm", plot = P4a)
  
}



# Size analysis -----------------------------------------------------------


#file path for the tsv file and create path for the tsv file
Size_path = '_fbqc_SizeMean.tsv'
Std_path = '_fbqc_SizeStd.tsv'


Size_list = sapply(experiment, function(x) paste0(file_path1,  x, '/fbqc/',  x, Size_path), USE.NAMES=FALSE)
file_name_list = sapply(experiment, function(x) paste0( x), USE.NAMES=FALSE)

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (e in experiment){
  f_list = paste0(file_path2, 'Size_plot/',  e, '_Size_plot/')
  dir_list = c(dir_list, f_list)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)


Size_table = list() 

for (i in c(1:length(Size_list))){
  db1 = read.table(file = Size_list[i], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db1$Exp = file_name_list[i]
  Size_table[[i]] <- db1
}

sdb = do.call(rbind, Size_table)

db1 =  sdb
db1$Cycle <- gsub('Inc', '', db1$Cycle)
db1$Cycle <- as.integer(as.character(db1$Cycle))
db1 = db1 %>% select(-Tile.1)

dbm = db1  %>% group_by(Exp, Base, Cycle)  %>% summarise(QC.Value = mean(QC.Value))
dbm$Tile = "Mean"
dbme = db1  %>% group_by(Exp, Base, Cycle)  %>% summarise(QC.Value = median(QC.Value))
dbme$Tile = "Median"
db2 =  bind_rows(db1, dbm, dbme)
db2$Tile <- factor(db2$Tile, levels=mixedsort(unique(db2$Tile)))

db2$Base[db2$Base == "A"] <- "G1"
db2$Base[db2$Base == "T"] <- "G2"
db2$Base[db2$Base == "G"] <- "R3"
db2$Base[db2$Base == "C"] <- "R4"
#db2 <- rename(db2, Size = QC.Value)
db2$Tile <- factor(db2$Tile, levels=mixedsort(unique(db2$Tile)))


for (i in c(1 :length(file_name_list))){
  
  #Wrap plot
  num.row = numrow(db2)
  
  P1 = ggplot(db2 %>% filter(Exp == file_name_list[i]), aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=2, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text( size = x.font)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size=theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Cluster Size"), 
         x = "Cycle", y = "Size") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~. , nrow = num.row)
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'Size', '_wrap_plot', '.png'),
         width = 40, height = 20, units = "cm", plot = P1)
  
  
  
  #Mean plot
  P2 = ggplot(db2 %>% filter(Exp == file_name_list[i], Tile == "Mean"), aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=2, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text( size = 20)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = 20)) +
    theme(strip.text = element_text(size = 20)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Mean Cluster Size"), 
         x = "Cycle", y = "Size") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~. , ncol = 1)
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'Mean_Size', '.png'),
         width = 40, height = 20, units = "cm", plot = P2)
}


Std_list = sapply(experiment, function(x) paste0(file_path1,  x, '/fbqc/',  x, Std_path), USE.NAMES=FALSE)
file_name_list = sapply(experiment, function(x) paste0( x), USE.NAMES=FALSE)


Std_table = list() 

for (i in c(1:length(Std_list))){
  db1 = read.table(file = Std_list[i], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db1$Exp = file_name_list[i]
  Std_table[[i]] <- db1
}

sdb = do.call(rbind, Std_table)

db1 =  sdb
db1$Cycle <- gsub('Inc', '', db1$Cycle)
db1$Cycle <- as.integer(as.character(db1$Cycle))
db1 = db1 %>% select(-Tile.1)

dbm = db1  %>% group_by(Exp, Base, Cycle)  %>% summarise(QC.Value = mean(QC.Value))
dbm$Tile = "Mean"
dbme = db1  %>% group_by(Exp, Base, Cycle)  %>% summarise(QC.Value = median(QC.Value))
dbme$Tile = "Median"
db2 =  bind_rows(db1, dbm, dbme)
db2$Tile <- factor(db2$Tile, levels=mixedsort(unique(db2$Tile)))

db2$Base[db2$Base == "A"] <- "G1"
db2$Base[db2$Base == "T"] <- "G2"
db2$Base[db2$Base == "G"] <- "R3"
db2$Base[db2$Base == "C"] <- "R4"
#db2 <- rename(db2, Size = QC.Value)
db2$Tile <- factor(db2$Tile, levels=mixedsort(unique(db2$Tile)))


for (i in c(1 :length(file_name_list))){
  
  #Wrap plot
  num.row = numrow(db2)
  
  P1 = ggplot(db2 %>% filter(Exp == file_name_list[i]), aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=2, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text( size = x.font)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size=theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Cluster Size Stdev"), 
         x = "Cycle", y = "Size") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~. , nrow = num.row)
  
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'Size_Stdev', '_wrap_plot', '.png'),
         width = 40, height = 20, units = "cm", plot = P1)
  
  
  
  #Mean plot
  P2 = ggplot(db2 %>% filter(Exp == file_name_list[i], Tile == "Mean"), aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=2, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text( size = 20)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = 20)) +
    theme(strip.text = element_text(size=20)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Mean Cluster Size Stdev"), 
         x = "Cycle", y = "Size") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile ~. , ncol = 1)
  
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'Mean_Size_Stdev', '.png'),
         width = 40, height = 20, units = "cm", plot = P2)
}


# Floursence Intensity Analysis -------------------------------------------

x.font = 13
y.font = 13
theme.font = 13
mx.font = 13


#file path for the tsv file and create path for the tsv file
relative_path = '_fbqc_RelativeIntensityMean.tsv'
background_path = '_fbqc_BackgroundMean.tsv'

relative_list = sapply(experiment, function(x) paste0(file_path1,  x, '/fbqc/',  x, relative_path), USE.NAMES=FALSE)
background_list = sapply(experiment, function(x) paste0(file_path1,  x, '/fbqc/',  x, background_path), USE.NAMES=FALSE)
file_name_list = sapply(experiment, function(x) paste0( x), USE.NAMES=FALSE)

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (e in experiment){
  f_list = paste0(file_path2, 'fluorescent_intensity_plot/',  e, '_fluorescent_intensity_plot/')
  dir_list = c(dir_list, f_list)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)


relative_table = list() 
background_table = list() 

for (i in c(1:length(relative_list))){
  db1 = read.table(file = relative_list[i], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db1 = db1 %>% filter(QC.Value != 0)
  db1$Exp = file_name_list[i]
  db1$FL.type = "Relative"
  relative_table[[i]] <- db1
}

rdb = do.call(rbind, relative_table)

for (i in c(1:length(background_list))){
  db2 = read.table(file = background_list[i], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db2 = db2 %>% filter(QC.Value != 0)
  db2$Exp = file_name_list[i]
  db2$FL.type = "Background"
  background_table[[i]] <- db2
}

bdb = do.call(rbind, background_table)

db1 =  bind_rows(rdb, bdb)
db1$Cycle <- gsub('Inc', '', db1$Cycle)
db1$Cycle <- as.integer(as.character(db1$Cycle))
db1 = db1 %>% select(-Tile.1)

dbm = db1  %>% group_by(Exp, Base, Cycle, FL.type)  %>% summarise(QC.Value = mean(QC.Value))
dbm$Tile = "Mean"
dbme = db1  %>% group_by(Exp, Base, Cycle, FL.type)  %>% summarise(QC.Value = median(QC.Value))
dbme$Tile = "Median"
db2 =  bind_rows(db1, dbm, dbme)
db2$Tile <- factor(db2$Tile, levels=mixedsort(unique(db2$Tile)))

db2$Base[db2$Base == "A"] <- "G1"
db2$Base[db2$Base == "T"] <- "G2"
db2$Base[db2$Base == "G"] <- "R3"
db2$Base[db2$Base == "C"] <- "R4"
colnames(db2)[2] = "FL"
db2$Tile_Type = paste(db2$Tile, db2$FL.type)

channel_list = c("G1", "G2", "R3", "R4")

datalist2 = list()


for (e in c(1 :length(file_name_list))) {
  for (i in c(1 :length(channel_list))) {
    r.vector = c()
    db3a = db2 %>% filter(Base == channel_list[i], Exp == file_name_list[e])
    db3 = db3a %>% group_by(Tile_Type)
    fitted_models =  db3 %>% do(model = lm(log(FL) ~ Cycle, data = .))
    co_db = fitted_models %>% ungroup() %>% transmute(Tile_Type, Coef = map(model, tidy)) %>% unnest(Coef) %>% 
      filter(term == "Cycle") %>% select(Tile_Type, estimate)
    for (c in c(1: nrow(fitted_models))) r.vector[c] <- summary(fitted_models$model[[c]])$r.squared
    co_db$r.squared = r.vector
    co_db$Base = channel_list[i]  
    datalist2[[paste0(channel_list[i],"_", experiment_list[e])]] <- co_db
  }
}


#Combind the exponential fit data 
#db4 = do.call(rbind, datalist2)
db4 = bind_rows( datalist2, .id = 'Exp')
db4$perchange = (exp(db4$estimate)-1)
db4 = db4 %>% separate(Tile_Type, c("Tile","FL.type"), sep = " ")
db4$Tile_Type = paste(db4$Tile, db4$FL.type)

#Create a column to specifify the experiment name
exp_col = gsub("G\\d_", "", db4$Exp)
exp_col = gsub("G\\d_", "", exp_col)
exp_col = gsub("R\\d_", "", exp_col)
exp_col = gsub("\\.\\d+", "", exp_col)
db4$Exp = exp_col
db4$Tile = gsub("mm", "", db4$Tile)
db4$Tile <- factor(db4$Tile, levels=mixedsort(unique(db4$Tile)))

for (i in c(1 :length(file_name_list))){
  
  dd = db2 %>% filter(Exp == file_name_list[i], Tile != "Median")
  tile.num  = length(unique(dd$Tile))
  
  if (tile.num <= 3){
    num.row = 1
  } else if (tile.num %in% c(4, 6)){
    num.row = 2
  } else if (tile.num %in% c(5)){
    num.row = 3
  } else if (tile.num %in% c(7, 8)){
    num.row = 4
  } else if (tile.num %in% c(9, 11, 10, 13, 14, 15)){
    num.row = 5
  } else if (tile.num %in% c( 16, 17, 18)){
    num.row = 6
  } 
  
  #Wrap plot
  P1 = ggplot(db2 %>% filter(Exp == file_name_list[i], Tile != "Median"), aes(x = Cycle, y = FL, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 1, alpha = 0.7) + 
    geom_point(size=1, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size= theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Fluorescent Intensity"), 
         x = "Cycle", y = "Fluorescent Intensity") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile_Type ~. , nrow = num.row)
  
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'Fluorescent_Intensity', '_wrap_plot', '.png'),
         width = 40, height = 20, units = "cm", plot = P1)
  
  P1a = ggplot(db2 %>% filter(Exp == file_name_list[i], Tile != "Median"), aes(x = Cycle, y = FL, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 1, alpha = 0.7) + 
    geom_point(size=1, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size= theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Fluorescent Intensity"), 
         x = "Cycle", y = "Fluorescent Intensity") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile_Type ~. , nrow = num.row, scales = "free_y")
  
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'Fluorescent_Intensity', '_wrap_plot_2', '.png'),
         width = 40, height = 20, units = "cm", plot = P1a)
  
  #Mean plot
  P2 = ggplot(db2 %>% filter(Exp == file_name_list[i], Tile == "Mean"), aes(x = Cycle, y = FL, group = Base, colour=Base)) + 
    theme_bw() +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=2, alpha = 0.7) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text( size = 20)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
          axis.text.y = element_text(hjust = 1, size = 18)) +
    theme(strip.text = element_text(size=14)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(file_name_list[i], " Fluorescent Intensity"), 
         x = "Cycle", y = "Fluorescent Intensity") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    facet_wrap(Tile_Type ~. , ncol = 2)
  
  
  ggsave(filename = paste0(dir_list[i],  experiment[i], '_', 'Mean_Fluorescent_Intensity', '.png'),
         width = 40, height = 20, units = "cm", plot = P2)
  
  
  for (i in c(1:length(file_name_list))){
    P3 = ggplot(db4 %>% filter(Exp == file_name_list[i]), aes(x = Tile, y = perchange, colour= Base, group = Base)) + 
      theme_bw() +
      geom_line(size = 1.2, alpha = 0.7) + 
      geom_point(size=4, alpha = 0.7) + 
      scale_y_continuous(labels=scales::percent) +
      labs(title = paste0(file_name_list[i], " Fluorescent Intensity intensity percentage change per cycle"), x = "Tile", y = "Percentage change per cycle") +
      scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6"),
                         name = "Filter") +
      theme(plot.title = element_text(hjust = 0.5, size=25)) +
      theme(axis.title.x = element_text(size=25), axis.text.x = element_text(angle = 45, hjust = 0.9, size = mx.font)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=25), 
            axis.text.y = element_text(size = 22)) +
      theme(legend.title = element_text(size = 25), legend.text = element_text(size = 22)) +
      theme(strip.text = element_text(size=25)) +
      theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
      facet_grid(.~ FL.type)
    
    
    ggsave(filename = paste0(dir_list[i], file_name_list[i], '_', "intensity_change_per_cycle", '.png'),
           width = 40, height = 20, units = "cm", plot = P3)
  }
}

#Analysis for only incorporation


# Phasing -----------------------------------------------------------------


experiment = experiment[str_detect(experiment, "Inc")]


#file path for the tsv file and create path for the tsv file
base_path = '_phase_base.tsv'
cycle_path = '_phase_cycle.tsv'
cross2_path = '_phase_and_crosstalk_2_cycle.tsv'

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (e in experiment){
  path1 = paste0(file_path1,  e)
  d_list = list.dirs(path1, recursive = FALSE, full.names = FALSE)
  p_list = subset(d_list, grepl("param_", d_list))
  path2 = 'C:/Users/kelvin/Documents/Figures/Phasing_plot/'
  f_list = sapply(p_list, function(x) paste0(path2,  e, '_phase_plot/', x, '/'), 
                  USE.NAMES=FALSE)
  dir_list = c(dir_list, f_list)
}

base_list = c()
for (e in experiment){
  path1 = paste0(file_path1,  e)
  d_list = list.dirs(path1, recursive = FALSE, full.names = FALSE)
  p_list = subset(d_list, grepl("param_", d_list))
  f_list = sapply(p_list, function(x) paste0(path1, '/', x, '/phase/',  e, base_path), 
                  USE.NAMES=FALSE)
  base_list = c(base_list, f_list)
}

cycle_list = c()
for (e in experiment){
  path1 = paste0(file_path1,  e)
  d_list = list.dirs(path1, recursive = FALSE, full.names = FALSE)
  p_list = subset(d_list, grepl("param_", d_list))
  f_list = sapply(p_list, function(x) paste0(path1, '/', x, '/phase/',  e, cycle_path), 
                  USE.NAMES=FALSE)
  cycle_list = c(cycle_list, f_list)
}

cross2_list = c()
for (e in experiment){
  path1 = paste0(file_path1,  e)
  d_list = list.dirs(path1, recursive = FALSE, full.names = FALSE)
  p_list = subset(d_list, grepl("param_", d_list))
  f_list = sapply(p_list, function(x) paste0(path1, '/', x, '/phase/',  e, cross2_path), 
                  USE.NAMES=FALSE)
  cross2_list = c(cross2_list, f_list)
}

file_name =c()
para_name = c()
for (f in base_list){
  p_name = unlist(rm_between(f, '_Inc/', '/phase', extract=TRUE))
  f_name = unlist(rm_between(f, 'Bioinformatics_analysis/', '/', extract=TRUE))
  p_name = gsub('CRT', '', p_name)
  file_name = c(file_name, f_name)
  para_name = c(para_name, p_name)
}

base_list = unlist(base_list)
cycle_list = unlist(cycle_list)
cross2_list = unlist(cross2_list)

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)

#Bases plot
#Plotting each graph, mean, median graph, wrap graph
for (f in c(1:length(base_list))){
  db = read.delim(file = base_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db$PP1 = sapply(db$PP, function(x) str_extract(x, "\\w?\\d\\d\\.?\\d?\\d?mm"))
  db$Tile = db$PP1
  num.row = numrow(db)
  
  for (i in unique(db$PP)){
    P1 = ggplot(db %>% filter(PP == i), aes(x = Cycle, y = Correlation, colour=Base, group = Base)) + 
      theme_bw() +
      geom_line(size = 2, alpha = 0.5) + 
      geom_point(size=3, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', ' ', i), x = "Cycle", y = "Correlation") +
      scale_x_continuous(breaks= pretty_breaks()) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 28)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], '_', i,'.png'), 
           width = 40, height = 20, units = "cm", plot = P1)
  }
  
  P2 = ggplot(db, aes(x = Cycle, y = Correlation, colour=Base, group = Base)) + 
    theme_bw() +
    geom_line(size = 0.7, alpha = 0.5) + 
    geom_point(size=0.5, alpha = 0.7) + 
    labs(title = paste0(file_name[f], ' (', para_name[f], ')', ' ', "Correlations Between Bases"), x = "Cycle", y = "Correlation") +
    scale_x_continuous(breaks= pretty_breaks()) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text( vjust=0.6, angle = 50, size = x.font))+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = y.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size= theme.font)) +
    theme(panel.spacing=unit(2, "lines")) +
    facet_wrap(PP1~., nrow = num.row)  
  
  ggsave(filename = paste0(dir_list[f], file_name[f], '_', "correlation_between_bases", '_wrap_plot', '.png'), 
         width = 40, height = 20, units = "cm", plot = P2)
  
  
  P3 = ggplot(db %>% group_by(Cycle, Base) %>% summarise(mean_corr = mean(Correlation)), 
              aes(x = Cycle, y = mean_corr, colour=Base, group = Base)) + 
    theme_bw() +
    geom_line(size = 2, alpha = 0.5) + 
    geom_point(size=3, alpha = 0.7) + 
    labs(title = paste0(file_name[f], ' (', para_name[f], ')', ' ', "Correlation between bases (Mean)"), x = "Cycle", y = "Correlation") +
    scale_x_continuous(breaks= pretty_breaks()) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 28)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = 28)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size=28)) +
    theme(plot.margin = margin(2, 2, 2, 2, "cm"))
  
  ggsave(filename = paste0(dir_list[f], file_name[f], '_', "Correlation_between_bases", '_mean_plot', '.png'), 
         width = 40, height = 20, units = "cm", plot = P3)
  
  P3b = ggplot(db %>% group_by(Cycle, Base) %>% summarise(median_corr = median(Correlation)), 
               aes(x = Cycle, y = median_corr, colour=Base, group = Base)) + 
    theme_bw() +
    geom_line(size = 2, alpha = 0.5) + 
    geom_point(size=3, alpha = 0.7) + 
    labs(title = paste0(file_name[f], ' (', para_name[f], ')', ' ', "Correlation between bases (Median)"), x = "Cycle", y = "Correlation") +
    scale_x_continuous(breaks= pretty_breaks()) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 28)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = 28)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size=28)) +
    theme(plot.margin = margin(2, 2, 2, 2, "cm"))
  
  ggsave(filename = paste0(dir_list[f], file_name[f], '_', "Correlation_between_bases", '_median_plot', '.png'), 
         width = 40, height = 20, units = "cm", plot = P3b)
  
}

#Cycle plot
for (f in c(1:length(cycle_list))){
  db = read.delim(file = cycle_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db$PP1 = sapply(db$PP, function(x) str_extract(x, "\\w?\\d\\d\\.?\\d?\\d?mm"))
  db$Cycle <- factor(db$Cycle, levels=mixedsort(unique(db$Cycle)))
  labels <- db$Cycle[ sort(c(seq(1, nrow(db)/length(unique(db$PP)), by= 20), length(unique(db$Cycle))*4 ))]
  labels2 <- db$Cycle[ sort(c(seq(1, nrow(db)/length(unique(db$PP)), by= 40), length(unique(db$Cycle))*4 ))]
  
  db$Tile = db$PP1
  num.row = numrow(db)
  
  for (i in unique(db$PP)){
    P4 = ggplot(db %>% filter(PP == i), aes(x = Cycle, y = Correlation, colour=Base, group = Base)) + 
      theme_bw() +
      geom_line(size = 2, alpha = 0.5) + 
      geom_point(size=3, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', ' ', i), x = "Cycle", y = "Correlation") +
      scale_x_discrete(breaks=labels) +
      scale_color_manual(values = c("A" = "#3cb44b", "T" = "#f58231", "G" = "#4363d8", "C" = "#f032e6")) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20, vjust=0.6, angle = 50)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], '_', i,'.png'), 
           width = 40, height = 20, units = "cm", plot = P4)
  }
  
  P5 = ggplot(db %>% filter(Cycle != "All"), aes(x = Cycle, y = Correlation, colour=Base, group = Base)) + 
    theme_bw() +
    geom_line(size = 0.7, alpha = 0.5) + 
    geom_point(size=0.5, alpha = 0.7) + 
    labs(title = paste0(file_name[f], ' (', para_name[f], ')', ' ', "Correlations Between Cycles"), x = "Cycle", y = "Correlation") +
    scale_x_discrete(breaks=labels2) +
    scale_color_manual(values = c("A" = "#3cb44b", "T" = "#f58231", "G" = "#4363d8", "C" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text( vjust=0.6, angle = 50, size = x.font))+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = y.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size= theme.font)) +
    theme(panel.spacing=unit(2, "lines")) +
    facet_wrap(PP1~., nrow = num.row)  
  
  ggsave(filename = paste0(dir_list[f], file_name[f], '_', "Correlation_between_cycles", '_wrap_plot', '.png'), 
         width = 40, height = 20, units = "cm", plot = P5)
  
  P6 = ggplot(db %>% group_by(Cycle, Base) %>% summarise(mean_corr = mean(Correlation)), 
              aes(x = Cycle, y = mean_corr, colour=Base, group = Base)) + 
    theme_bw() +
    geom_line(size = 2, alpha = 0.5) + 
    geom_point(size=3, alpha = 0.7) + 
    labs(title = paste0(file_name[f], ' (', para_name[f], ')', ' ', "Correlation between cycles (Mean)"), x = "Cycle", y = "Correlation") +
    scale_x_discrete(breaks=labels) +
    scale_color_manual(values = c("A" = "#3cb44b", "T" = "#f58231", "G" = "#4363d8", "C" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20, vjust=0.6, angle = 50)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = 28)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size=28)) +
    theme(plot.margin = margin(2, 2, 2, 2, "cm"))
  
  ggsave(filename = paste0(dir_list[f], file_name[f], '_', "Correlation_between_cycles", '_mean_plot', '.png'), 
         width = 40, height = 20, units = "cm", plot = P6)
  
  
  P6b = ggplot(db %>% group_by(Cycle, Base) %>% summarise(median_corr = median(Correlation)), 
               aes(x = Cycle, y = median_corr, colour=Base, group = Base)) + 
    theme_bw() +
    geom_line(size = 2, alpha = 0.5) + 
    geom_point(size=3, alpha = 0.7) + 
    labs(title = paste0(file_name[f], ' (', para_name[f], ')', ' ', "Correlation between cycles (Median)"), x = "Cycle", y = "Correlation") +
    scale_x_discrete(breaks=labels) +
    scale_color_manual(values = c("A" = "#3cb44b", "T" = "#f58231", "G" = "#4363d8", "C" = "#f032e6")) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20, vjust=0.6, angle = 50)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = 28)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size=28)) +
    theme(plot.margin = margin(2, 2, 2, 2, "cm"))
  
  ggsave(filename = paste0(dir_list[f], file_name[f], '_', "Correlation_between_cycles", '_median_plot', '.png'), 
         width = 40, height = 20, units = "cm", plot = P6b)
  
}

#Cross_plot
for (f in c(1:length(cross2_list))){
  db = read.delim(file = cross2_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db$PP1 = sapply(db$PP, function(x) str_extract(x, "\\w\\d\\d\\.\\d\\dmm"))
  db$Cycles <- factor(db$Cycles, levels=mixedsort(unique(db$Cycles)))
  labels <- db$Cycles[ sort(c(seq(1, nrow(db)/length(unique(db$PP)), by= 10), length(unique(db$Cycles))*2 ))]
  labels2 <- db$Cycles[ sort(c(seq(1, nrow(db)/length(unique(db$PP)), by= 20), length(unique(db$Cycles))*2 ))]
  
  db$Tile = db$PP1
  num.row = numrow(db)
  
  for (i in unique(db$PP)){
    P1 = ggplot(db %>% filter(PP == i), aes(x = Cycles, y = Percentage, colour=QC, group = QC)) + 
      theme_bw() +
      geom_line(size = 2, alpha = 0.5) + 
      geom_point(size=3, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', ' ', i), x = "Cycle", y = "Percentage") +
      scale_x_discrete(breaks=labels) +
      scale_y_continuous(labels=scales::percent) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20, vjust=0.6, angle = 50)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], '_cross2_', i,'.png'), 
           width = 40, height = 20, units = "cm", plot = P1)
  }
  
  
  P2 = ggplot(db %>% filter(Cycles != "All"), aes(x = Cycles, y = Percentage, colour=QC, group = QC)) +
    theme_bw() +
    geom_line(size = 0.7, alpha = 0.5) + 
    geom_point(size=0.5, alpha = 0.7) + 
    labs(title = paste0(file_name[f], ' (', para_name[f], ')'), x = "Cycle", y = "Percentage") +
    scale_x_discrete(breaks=labels2) +
    scale_y_continuous(labels=scales::percent) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text( vjust=0.6, angle = 50, size = x.font))+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = y.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size= theme.font)) +
    theme(panel.spacing=unit(2, "lines")) +
    facet_wrap(PP1~., nrow = num.row)  
  
  
  ggsave(filename = paste0(dir_list[f], file_name[f], '_cross2_', '_wrap_plot','.png'), 
         width = 40, height = 20, units = "cm", plot = P2)
  
  
  P3 = ggplot(db %>% group_by(Cycles, QC) %>% summarise(mean_percent = mean(Percentage)), 
              aes(x = Cycles, y = mean_percent, colour=QC, group = QC)) + 
    theme_bw() +
    geom_line(size = 2, alpha = 0.5) + 
    geom_point(size=3, alpha = 0.7) + 
    labs(title = paste0(file_name[f], ' (', para_name[f], ')', " Phasing and Cross Talking (Mean)"), x = "Cycles", y = "Percentage") +
    scale_x_discrete(breaks=labels) +
    scale_y_continuous(labels=scales::percent) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20, vjust=0.6, angle = 50)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = 28)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size=28)) +
    theme(plot.margin = margin(2, 2, 2, 2, "cm"))
  
  ggsave(filename = paste0(dir_list[f], file_name[f],  "_cross2_", '_mean_plot', '.png'), 
         width = 40, height = 20, units = "cm", plot = P3)
  
  
  
  P3b = ggplot(db %>% group_by(Cycles, QC) %>% summarise(median_percent = median(Percentage)), 
               aes(x = Cycles, y = median_percent, colour=QC, group = QC)) + 
    theme_bw() +
    geom_line(size = 2, alpha = 0.5) + 
    geom_point(size=3, alpha = 0.7) + 
    labs(title = paste0(file_name[f], ' (', para_name[f], ')', " Phasing and Cross Talking (Median)"), x = "Cycles", y = "Percentage") +
    scale_x_discrete(breaks=labels) +
    scale_y_continuous(labels=scales::percent) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20, vjust=0.6, angle = 50)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = 28)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size=28)) +
    theme(plot.margin = margin(2, 2, 2, 2, "cm"))
  
  ggsave(filename = paste0(dir_list[f], file_name[f],  "_cross2_", '_median_plot', '.png'), 
         width = 40, height = 20, units = "cm", plot = P3b)
}


for (f in c(1:length(cross2_list))){
  db = read.delim(file = cross2_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db1 = db %>% filter(Cycles != "All")
  db1$Cycles = gsub("-\\d+", "", db1$Cycles)
  db1$Cycles = as.numeric(db1$Cycles)
  db1$PP1 = sapply(db1$PP, function(x) str_extract(x, "\\w?\\d\\d\\.?\\d?\\d?mm"))
  
  db1$Tile = db1$PP1
  num.row = numrow(db1)
  
  for (i in unique(db$PP)){
    P4 = ggplot(db1 %>% filter(PP == i), aes(x = Cycles, y = Percentage, colour=QC, group = QC)) + 
      theme_bw() +
      geom_line(size = 1, alpha = 0.3) + 
      geom_point(size=3, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', ' ', i), x = "Cycle", y = "Percentage") +
      geom_smooth(method = "lm", se=FALSE, formula = y ~ x,  alpha = 0.3) +
      stat_poly_eq(formula = y ~ x, 
                   aes(label = paste(..eq.label.., sep = "~~~")), 
                   label.x = "left", label.y = "top", parse = TRUE, size = 10) +
      scale_x_continuous(breaks= pretty_breaks()) +
      scale_y_continuous(labels=scales::percent) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], '_cross2_linear_fit ', i,'.png'), 
           width = 40, height = 20, units = "cm", plot = P4)
  }
  
  P5 = ggplot(db1 %>% group_by(Cycles, QC) %>% summarise(mean_percent = mean(Percentage)), 
              aes(x = Cycles, y = mean_percent, colour=QC, group = QC)) + 
    theme_bw() +
    geom_line(size = 1, alpha = 0.3) + 
    geom_point(size=3, alpha = 0.7) + 
    geom_smooth(method = "lm", se=FALSE, formula = y ~ x,  alpha = 0.3) +
    stat_poly_eq(formula = y ~ x, 
                 aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x = "left", label.y = "top", parse = TRUE, size = 10) +
    labs(title = paste0(file_name[f], ' (', para_name[f], ')', " linear fit (mean)"), x = "Cycles", y = "Percentage") +
    scale_x_continuous(breaks= pretty_breaks()) +
    scale_y_continuous(labels=scales::percent) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = 28)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size=28)) +
    theme(plot.margin = margin(2, 2, 2, 2, "cm"))
  
  ggsave(filename = paste0(dir_list[f], file_name[f],  "_cross2_", '_linear_fit_mean_plot', '.png'), 
         width = 40, height = 20, units = "cm", plot = P5)
  
  P6 = ggplot(db1, aes(x = Cycles, y = Percentage, colour=QC, group = QC)) +
    theme_bw() +
    geom_line(size = 0.7, alpha = 0.5) + 
    geom_point(size=0.5, alpha = 0.7) + 
    geom_smooth(method = "lm", se=FALSE, formula = y ~ x,  alpha = 0.3) +
    stat_poly_eq(formula = y ~ x, 
                 aes(label = paste(..eq.label.., sep = "~~~")), 
                 label.x.npc = "left",  parse = TRUE, size = 5) +
    labs(title = file_name[f], ' (', para_name[f], ')', " linear fit (mean)", x = "Cycles", y = "Percentage") +
    scale_x_continuous(breaks= pretty_breaks()) +
    scale_y_continuous(labels=scales::percent) +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=30), axis.text.x = element_text( vjust=0.6, angle = 50, size = x.font))+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
          axis.text.y = element_text(size = y.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    theme(strip.text = element_text(size= theme.font)) +
    theme(panel.spacing=unit(2, "lines")) +
    facet_wrap(PP1~., nrow= num.row)  
  
  
  ggsave(filename = paste0(dir_list[f], file_name[f],  "_cross2_", '_linear_fit_mean_plot_wrap', '.png'), 
         width = 40, height = 20, units = "cm", plot = P6)
  
}  

#Phase_percent
#file path for the tsv file and create path for the tsv file

cross_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = '_phase_and_crosstalk_\\d+_cycle.tsv$', full.names = TRUE, recursive = TRUE)
  cross_list = c(cross_list, x)
}

phase.percent.list = cross_list[!grepl('_Inc_phase_and_crosstalk_2', cross_list)]
phase.percent.list2 = cross_list[grepl('_Inc_phase_and_crosstalk_2', cross_list)]
para_name = unlist(rm_between(phase.percent.list, '_Inc/', '/phase', extract=TRUE))
file_name = unlist(rm_between(phase.percent.list, 'Bioinformatics_analysis/', '/', extract=TRUE))

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (i in c(1: length(file_name))){
  x = paste0('C:/Users/kelvin/Documents/Figures/Phasing_plot/', file_name[i], '_phase_plot/', para_name[i], '/')
  dir_list = c(dir_list, x)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)

#
for (f in c(1:length(phase.percent.list))){
  df = read.delim(file = phase.percent.list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  df2 = read.delim(file = phase.percent.list2[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  df = df[complete.cases(df), ]
  df2 = df2[complete.cases(df2), ]
  df = filter(df, Percentage != 0)
  df2 = filter(df2, Percentage != 0)
  db = rbind(df, df2)
  db$PP1 = sapply(db$PP, function(x) str_extract(x, "\\w?\\d\\d\\.?\\d?\\d?mm"))
  db = db %>% filter(Cycles != "All")
  db = db %>% separate(Cycles, c("Cycles", "Cycles2"))
  db$QC = gsub("Phasing\\+", "Phasing", db$QC)
  db$QC = gsub("Phasing-", "Prephasing", db$QC)
  db$PP1 = gsub("mm", "", db$PP1)
  db$Cycles = as.numeric(db$Cycles)
  db1 = db %>% filter(Cycles <= 25, QC %in% c("Phasing Off-Diag", "Prephasing Off-Diag", "Phasing Prob", "Prephasing Prob")) %>% 
    group_by(QC, PP1) %>% summarise(mean_phase = mean(Percentage)) 
  db2 = db1 %>% summarise(mean_phase = mean(mean_phase))
  db2$ PP1 = "mean"
  db3 = bind_rows(db1, db2)
  db3$QC = gsub("Diag", "Diag (mean % of 25 cycles)", db3$QC)
  db3$QC = gsub("Prob", "Prob (mean % of 25 cycles)", db3$QC)
  
  db4 = db %>% filter(QC %in% c("Phasing Prob", "Prephasing Prob")) %>% 
    group_by(QC, PP1) %>% summarise(mean_phase = mean(Percentage)) 
  db5 = db4 %>% summarise(mean_phase = mean(mean_phase))
  db5$PP1 = "mean"
  db6 = bind_rows(db4, db5)
  db6$QC = gsub("Prob", "Prob (mean % of All cycles)", db6$QC)
  
  db7 = bind_rows(db3, db6)
  
  
  db8 = db %>% group_by(QC, Cycles) %>% summarise(Percentage = mean(Percentage))
  db8$PP1 = "mean"
  db9 = db %>% select(QC, Cycles, Percentage, PP1)
  db10 = bind_rows(db8, db9)
  db10$Cycles = as.numeric(db10$Cycles)
  
  P7 = ggplot(db10 %>% filter(QC %in% c("Phasing Off-Diag", "Prephasing Off-Diag", "Phasing Prob", "Prephasing Prob")),
              aes(x = Cycles, y =Percentage, col= QC, group = QC)) + 
    theme_bw() +
    geom_line(size = 1, alpha = 0.7) + 
    geom_point(size=1.5, alpha = 0.7) + 
    scale_y_continuous(labels=scales::percent) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = "Cycle", y = "Percentage",
         title = paste0(file_name[f], ' (', para_name[f], ')')) +
    scale_color_manual(values = c("Prephasing Off-Diag" = "orange", 
                                  "Phasing Off-Diag" = "red", 
                                  "Prephasing Prob" = "green", 
                                  "Phasing Prob" = "blue")) +
    theme(plot.title = element_text(hjust = 0.5, size=20)) +
    theme(axis.title.x=element_text(size=20), axis.text.x = element_text(size = 13)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=20),
          axis.text.y = element_text(size = 13)) +
    theme(legend.title = element_text(size = 20), legend.text=element_text(size=15)) +
    theme(strip.text = element_text(size=15)) +
    theme(plot.margin = margin(0.25, 1, 1, 1, "cm")) +
    theme(legend.position="bottom") +
    guides(col = guide_legend(nrow = 1)) +
    facet_wrap(PP1 ~., nrow = 2 )
  
  if (length(unique(db10$PP1)) < 15) {
    f.size = 15} else if (length(unique(db10$PP1)) >= 15) {
      f.size = 14
    }
  
  db7 = db7 %>% spread(PP1, mean_phase)
  db7 = db7[c(6,3,5,2, 4,1),]
  db7[, -1] = db7[, -1]*100
  db7 <- db7 %>% 
    mutate_if(is.numeric, round, digits = 2)
  
  tt3 <- ttheme_minimal(
    core=list(bg_params = list(blues9[1:6], col=NA),
              fg_params=list(fontface=2, fontsize = f.size)),
    colhead=list(fg_params=list(col="navyblue", fontface=4, fontsize = f.size)))
  
  tbl =  tableGrob(db7, rows =  NULL, theme = tt3)
  
  P7a = grid.arrange(tbl, P7,  
                     nrow=2, 
                     heights=c(3, 8))
  
  ggsave(filename = paste0(dir_list[f], file_name[f], '_', para_name[f], '_phase_percent','.png'), 
         width = 40, height = 20, units = "cm", plot = P7a)
}  


# Summary report ----------------------------------------------------------


summary_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = '_summary_report.tsv$', full.names = TRUE, recursive = TRUE)
  summary_list = c(summary_list, x)
}


para_name = unlist(rm_between(summary_list, '_Inc/', '/analyze', extract=TRUE))
file_name = unlist(rm_between(summary_list, 'Bioinformatics_analysis/', '/', extract=TRUE))

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (i in c(1: length(file_name))){
  x = paste0('C:/Users/kelvin/Documents/Figures/Phasing_plot/', file_name[i], '_phase_plot/', para_name[i], '/')
  dir_list = c(dir_list, x)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)


for (f in c(1:length(summary_list))){
  df = read.delim(file = summary_list[f], sep = '\t', header = T, stringsAsFactors = FALSE)
  #df = df[c(-21), ]
  colnames(df) = gsub("X", "", colnames(df))
  colnames(df) = gsub("Category", "Run Statistics", colnames(df))
  df = df[c(-1),]
  
  a =  which(df$`Run Statistics` == "Read quality") -1
  b = a +1
  c = which(df$`Run Statistics` == "Alignment") - 1
  d = which(df$`Run Statistics` == "Alignment")
  
  df1 = df[ c(1:a), ]
  
  df2 = df[ c(b:c), ]
  colnames(df2) = gsub("Run Statistics", "Read quality", colnames(df2))
  df2 = df2[c(-1),]
  
  df3 = df[ c(d:nrow(df)), ]
  colnames(df3) = gsub("Run Statistics", "Alignment", colnames(df3))
  df3 = df3[c(-1),]
  
  if (ncol(df) == 2){
    base.size = 24
  } else if (ncol(df) == 3){
    base.size = 22
  } else if (ncol(df) == 4){
    base.size = 20
  } else if (ncol(df) >= 5){
    base.size = 16
  } 
  
  tt <- ttheme_minimal(base_size = base.size, colhead=list(fg_params=list(col="navyblue", fontface=4L)))
  
  tb1 = tableGrob(df1, rows =  NULL, theme=tt)
  tb1 <- gtable_add_grob(tb1,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 2, b = nrow(tb1), l = 1, r = ncol(tb1))
  tb1 <- gtable_add_grob(tb1,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, l = 1, r = ncol(tb1))
  
  tb2 = tableGrob(df2, rows =  NULL, theme=tt)
  tb2 <- gtable_add_grob(tb2,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 2, b = nrow(tb2), l = 1, r = ncol(tb2))
  tb2 <- gtable_add_grob(tb2,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, l = 1, r = ncol(tb2))
  
  tb3 = tableGrob(df3, rows =  NULL, theme=tt)
  tb3 <- gtable_add_grob(tb3,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 2, b = nrow(tb3), l = 1, r = ncol(tb3))
  tb3 <- gtable_add_grob(tb3,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, l = 1, r = ncol(tb3))
  
  haligned <- gtable_combine(tb1, tb2, tb3, along=1)
  
  p = grid.arrange(haligned, nrow=1)
  
  ggsave(filename = paste0(dir_list[f], file_name[f], '_', para_name[f], '_summary_report','.png'),
         width = 60, height = 30, units = "cm", plot = p)
}



