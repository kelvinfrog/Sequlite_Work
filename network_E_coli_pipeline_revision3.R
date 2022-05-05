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
library(uniftest)
library(ggpubr)
library(Rsamtools)


# Input -------------------------------------------------------------------

#experiment number, only the number with NO CRT, but it needs to have b or t if the file name has it at
#the beginning
experiment  = c("CRT75y_Inc")
#experiment  = c("b34x_Inc", "CRT146x_Inc", "147x_Inc")

#file_path1 = '//sequlitenas/SequLite Storage/Chemistry/Experiment_Data/'
file_path1 = '//sequlitenas/SequLite Storage/Bioinformatics/'
file_path2 = '//sequlitenas/SequLite Storage/Chemistry/'
file_path_coverage = '//sequlitenas/SequLite Storage/Bioinformatics/'
  
exp = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1), pattern = e, full.names = TRUE, recursive = F)
  exp = c(exp, x)
}

experiment = gsub(file_path1, "", exp)

x.font = 10
y.font = 10
theme.font = 10
legend.row = 2

# customize the fonts for tiles numbers
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
  } else if (x > 20 & x <= 30) {
    row = 5
  } else if (x > 30) {
    row = 8 
  }
  return(row)
}

tablefont = function(df){
  if (length(unique(df$Tile)) < 15) {
    f.size = 15} else if (length(unique(df$Tile)) %in% c(15, 16)) {
      f.size = 13} else if (length(unique(df$Tile))%in% c(17, 18)) {
        f.size = 12} else if (length(unique(df$Tile)) %in% c(19, 20)) {
          f.size = 11} else if (length(unique(df$Tile)) > 19) {
            f.size = 10}
  return(f.size)
}

positionname  = function(letter){
  if (letter == "b"){
    letter = "_bottom_"
  } else if (letter == "t") {
    letter = "_top_"
  }
  return(letter)
}


# SNR, FOG, Size, Pixel ----------------------------------------------------------

type_path = c('_fbqc_SNRMean.tsv', '_fbqc_Fractionofgoodfeatures.tsv', 
              '_fbqc_SizeMean.tsv', '_fbqc_SizeStd.tsv', '_fbqc_PixelCountMean.tsv')
plot_path = c('SNR_plot/', 'FOG_plot/', 'Size_plot/', 'Size_plot/', 'Pixel_plot/')
plotname = c('SNR', 'FOG', 'Size', 'Standard_Deviation', 'Pixel')
plot_title = c('Signal to Noise Ratio', 'Fraction of Good Features', 'Cluster Size', 
               'Cluster Size Standard Deviation', 'Pixel Intensity (75 percentile)')

for (num in c(1 : length(type_path))){
  filelist = c()
  for (e in experiment) {
    x = list.files(path =  paste0(file_path1,  e),
                   pattern = type_path[num], full.names = TRUE, recursive = TRUE)
    filelist = c(filelist, x)
  }
  
  file_name = unlist(rm_between(filelist, file_path1, '/', extract=TRUE))
  
  #creat path for plot directory, tsv directory, filename 
  dir_list = c()
  for (i in c(1: length(file_name))){
    x = paste0(file_path2, plot_path[num], file_name[i], '_', plot_path[num])
    dir_list = c(dir_list, x)
  }
  
  #Create directory
  sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), 
         USE.NAMES=FALSE)
  
  for (f in c(1 : length(filelist))){
    
    df = read.delim(file = filelist[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>%
      drop_na() %>% filter(Tile != "Mean")
    
    for (p in unique(str_extract(df$Tile, "^\\w"))){
      df1 = df  %>% filter(str_detect(Tile, p)) 
      dfm = df1 %>%  group_by(Cycle, Base) %>% summarise(QC.Value = mean(QC.Value)) %>%
        mutate(Tile = "mean")
      dfmed = df1 %>%  group_by(Cycle, Base) %>% summarise(QC.Value = median(QC.Value)) %>%
        mutate(Tile = "median")
      df2 = bind_rows(df1, dfm, dfmed) %>%
        mutate(Base=recode(Base, "A" = "G1", "T" = "G2", "G" = "R3", "C" = "R4")) 
      
      num.row = numrow(df2)
      
      P1 = ggplot(df2, aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
        theme_bw() +
        geom_line(size = 0.9, alpha = 0.7) + 
        geom_point(size=0.5, alpha = 0.7) + 
        scale_y_continuous(breaks = scales::pretty_breaks()) +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231",
                                      "R3" = "#4363d8", "R4" = "#f032e6")) +
        theme(plot.title = element_text(hjust = 0.5, size=30)) +
        theme(axis.title.x=element_text(size=28), 
              axis.text.x = element_text(size = x.font)) +
        theme(axis.title.y = element_text(size = 28), 
              axis.text.y = element_text(hjust = 1, size = y.font)) +
        theme(strip.text = element_text(size= theme.font)) +
        theme(legend.title = element_text(size = 30), 
              legend.text=element_text(size=28)) +
        labs(title = paste0(file_name[f], positionname(p), plot_title[num]), 
             x = "Cycle", y = plotname[num], color = "Filter") +
        theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
        facet_wrap(Tile ~. , nrow = num.row)
      
      ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), plotname[num], 
                               '_wrap_plot', '.png'), 
             width = 40, height = 20, units = "cm", plot = P1)
      
      
      #Mean plot
      P2 = ggplot(df2 %>% filter(Tile == "mean"), 
                  aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
        theme_bw() +
        geom_line(size = 1.2, alpha = 0.7) + 
        geom_point(size=2, alpha = 0.7) + 
        scale_y_continuous(breaks = scales::pretty_breaks()) +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231",
                                      "R3" = "#4363d8", "R4" = "#f032e6")) +
        theme(plot.title = element_text(hjust = 0.5, size=30)) +
        theme(axis.title.x=element_text(size=28), 
              axis.text.x = element_text( size = 20)) +
        theme(axis.title.y = element_text(size = 28), 
              axis.text.y = element_text(hjust = 1, size = 18)) +
        theme(strip.text = element_text(size=20)) +
        theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
        labs(title = paste0(file_name[f], positionname(p), "Mean ", plot_title[num]), 
             x = "Cycle", y = plotname[num], color = "Filter") +
        theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 
      
      ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), plotname[num], 
                               '_Mean', '.png'), 
             width = 40, height = 20, units = "cm", plot = P2)
      
      if (type_path[num] == '_fbqc_PixelCountMean.tsv'){
        r.vector = c()
        fitted_models = df2 %>% 
          group_by(Tile, Base) %>% 
          do(model = lm(log(QC.Value) ~ Cycle, data = .)) 
        for (c in c(1: nrow(fitted_models))) r.vector[c] <- summary(fitted_models$model[[c]])$r.squared
        co_db = fitted_models %>% ungroup() %>% transmute(Tile, Coef = map(model, tidy)) %>% 
          unnest(Coef) %>% 
          filter(term == "Cycle") %>% select(Tile, estimate) %>% 
          mutate(r.squared = r.vector) %>% 
          mutate(Base = fitted_models$Base) %>% 
          mutate(perchange = exp(estimate)-1) %>% 
          mutate(Tile = factor(Tile, levels=mixedsort(unique(Tile))))
        
        P3 = ggplot(co_db, aes(x = Tile, y = perchange, color=Base, group = Base)) + 
          theme_bw() +
          geom_line(size = 1.2, alpha = 0.7) + 
          geom_point(size=4, alpha = 0.7) + 
          scale_y_continuous(labels=scales::percent) +
          labs(title = paste0(file_name[f], positionname(p), "Pixel intensity percentage change per cycle"),
               x = "Tile", y = "Percentage change per cycle", color = "Filter") +
          scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", 
                                        "R3" = "#4363d8", "R4" = "#f032e6"))+
          theme(plot.title = element_text(hjust = 0.5, size=25)) +
          theme(axis.title.x=element_text(size=25), 
                axis.text.x = element_text(size = 18, angle = 45, hjust = 0.9)) +
          theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),size=25), 
                axis.text.y = element_text(size = 22)) +
          theme(legend.title = element_text(size = 25), legend.text = element_text(size = 22)) +
          theme(strip.text = element_text(size=25))
        
        db7 = co_db %>% select(Base, perchange, Tile)
        names(db7) = c("Filter", "Percent_change", "Tile")
        db7 = db7 %>% spread(Tile, Percent_change)
        db7[, -1] = db7[, -1]*100
        db7 <- db7 %>% 
          mutate_if(is.numeric, round, digits = 2)
        
        f.size = tablefont(co_db)
        
        tt3 <- ttheme_minimal(
          core=list(bg_params = list(fill = c("#3cb44b", "#f58231", "#4363d8", "#f032e6"), col=NA),
                    fg_params=list(fontface=2, fontsize = f.size)),
          colhead=list(fg_params=list(col="navyblue", fontface=4, fontsize = f.size)))
        
        tbl =  tableGrob(db7, rows =  NULL, theme = tt3)
        
        P3a = grid.arrange(P3, tbl, 
                           nrow=2, 
                           heights=c(8, 2))
        
        ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), 'Pixel_Percent_Change', '.png'),
               width = 40, height = 20, units = "cm", plot = P3a)
      }
    }
  }
}

# Floursence Intensity Analysis -------------------------------------------
x.font = 10
y.font = 10
theme.font = 12
mx.font = 10


#file path for the tsv file and create path for the tsv file
relative_path = '_fbqc_RelativeIntensityMean.tsv'
background_path = '_fbqc_BackgroundMean.tsv'


relative_list = c()
background_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = relative_path, full.names = TRUE, recursive = TRUE)
  y = list.files(path =  paste0(file_path1,  e),
                 pattern = background_path, full.names = TRUE, recursive = TRUE)
  relative_list = c(relative_list, x)
  background_list = c(background_list, y)
}

file_name = unlist(rm_between(relative_list, file_path1, '/', extract=TRUE))

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (e in experiment){
  f_list = paste0(file_path2, 'fluorescent_intensity_plot/',  e, '_fluorescent_intensity_plot/')
  dir_list = c(dir_list, f_list)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)


for (f in c(1 : length(relative_list))){
  
  dfr = read.delim(file = relative_list[f], sep = '\t', header = TRUE, 
                   stringsAsFactors = FALSE) %>% drop_na() %>% filter(Tile != "Mean")
  dfb = read.delim(file = background_list[f], sep = '\t', header = TRUE, 
                   stringsAsFactors = FALSE) %>% drop_na() %>% filter(Tile != "Mean")
  
  for (p in unique(str_extract(dfr$Tile, "^\\w"))){
    
    dfr1 = dfr  %>% filter(str_detect(Tile, p)) %>% mutate(QC = 'Relative')
    dfb1 = dfb %>% filter(str_detect(Tile, p)) %>% mutate(QC = 'Background')
    
    df = bind_rows(dfr1, dfb1)
    dfm = df %>%  group_by(Cycle, Base, QC) %>% 
      summarise(QC.Value = mean(QC.Value)) %>% 
      mutate(Tile = "mean")
    dfmed = df %>%  group_by(Cycle, Base, QC) %>% 
      summarise(QC.Value = median(QC.Value)) %>% 
      mutate(Tile = "median")
    df1 = bind_rows(df, dfm, dfmed) %>% 
      mutate(Base = recode(Base, "A" = "G1", "T" = "G2", "G" = "R3", "C" = "R4")) %>%
      mutate(label = paste0(Tile, " ", QC))
    
    r.vector = c()
    fitted_models = df1 %>% 
      group_by(Tile, Base, QC) %>% 
      do(model = lm(log(QC.Value) ~ Cycle, data = .)) 
    for (c in c(1: nrow(fitted_models))) r.vector[c] <- summary(fitted_models$model[[c]])$r.squared
    co_db = fitted_models %>% ungroup() %>% transmute(Tile, Coef = map(model, tidy)) %>% 
      unnest(Coef) %>% 
      filter(term == "Cycle") %>% select(Tile, estimate) %>% 
      mutate(r.squared = r.vector) %>% 
      mutate(Base = fitted_models$Base) %>%
      mutate(QC = fitted_models$QC) %>% 
      mutate(perchange = exp(estimate)-1) %>% 
      mutate(Tile = factor(Tile, levels=mixedsort(unique(Tile))))
    
    tile.num  = length(unique(df1$Tile))
    
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
    } else if (tile.num > 18) {
      num.row = 5
    }
    
    #Wrap plot
    P1 = ggplot(df1, aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
      theme_bw() +
      geom_line(size = 1, alpha = 0.7) + 
      geom_point(size=1, alpha = 0.7) + 
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", 
                                    "R3" = "#4363d8", "R4" = "#f032e6")) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
            axis.text.y = element_text(hjust = 1, size = y.font)) +
      theme(strip.text = element_text(size= theme.font)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      labs(title = paste0(file_name[f], positionname(p), " Fluorescent Intensity"), 
           x = "Cycle", y = "Fluorescent Intensity", color = "Filter") +
      theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
      facet_wrap(label ~. , nrow = num.row)
    
    ggsave(filename = paste0(dir_list[f],  file_name[f], positionname(p), 'Fluorescent_Intensity',
                             '_wrap_plot', '.png'),
           width = 40, height = 20, units = "cm", plot = P1)
    
    P2 = ggplot(df1, aes(x = Cycle, y = QC.Value, group = Base, colour=Base)) + 
      theme_bw() +
      geom_line(size = 1, alpha = 0.7) + 
      geom_point(size=1, alpha = 0.7) + 
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", 
                                    "R3" = "#4363d8", "R4" = "#f032e6")) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
            axis.text.y = element_text(hjust = 1, size = y.font)) +
      theme(strip.text = element_text(size= theme.font)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      labs(title = paste0(file_name[f], positionname(p), " Fluorescent Intensity"), 
           x = "Cycle", y = "Fluorescent Intensity", color = "Filter") +
      theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
      facet_wrap(label ~. , nrow = num.row, scale = "free_y")
    
    ggsave(filename = paste0(dir_list[f],  file_name[f], positionname(p), 'Fluorescent_Intensity',
                             '_wrap_plot_2', '.png'),
           width = 40, height = 20, units = "cm", plot = P2)
    
    #Mean plot
    P3 = ggplot(df1 %>% filter(Tile == "mean"), aes(x = Cycle, y = QC.Value, 
                                                    group = Base, color = Base)) + 
      theme_bw() +
      geom_line(size = 1.2, alpha = 0.7) + 
      geom_point(size=2, alpha = 0.7) + 
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", 
                                    "R3" = "#4363d8", "R4" = "#f032e6")) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=28), axis.text.x = element_text( size = 20)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 28), 
            axis.text.y = element_text(hjust = 1, size = 18)) +
      theme(strip.text = element_text(size=14)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      labs(title = paste0(file_name[f], positionname(p), "Mean Fluorescent Intensity"), 
           x = "Cycle", y = "Fluorescent Intensity") +
      theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
      facet_wrap(QC ~. , nrow = 1)
    
    ggsave(filename = paste0(dir_list[f],  file_name[f], positionname(p), 
                             'Mean_Fluorescent_Intensity', '.png'),
           width = 40, height = 20, units = "cm", plot = P3)
    
    P4 = ggplot(co_db, aes(x = Tile, y = perchange, color=Base, group = Base)) + 
      theme_bw() +
      geom_line(size = 1.2, alpha = 0.7) + 
      geom_point(size=4, alpha = 0.7) + 
      scale_y_continuous(labels=scales::percent) +
      labs(title = paste0(file_name[f], positionname(p), 
                          'Fluorescent Intensity intensity percentage change per cycle'),
           x = "Tile", y = "Percentage change per cycle", color = "Filter") +
      scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", 
                                    "R3" = "#4363d8", "R4" = "#f032e6"))+
      theme(plot.title = element_text(hjust = 0.5, size=25)) +
      theme(axis.title.x=element_text(size=25), 
            axis.text.x = element_text(size = 18, angle = 45, hjust = 0.9)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),size=25), 
            axis.text.y = element_text(size = 22)) +
      theme(legend.title = element_text(size = 25), legend.text = element_text(size = 22)) +
      theme(strip.text = element_text(size=25)) +
      facet_wrap(QC ~. , nrow = 1)
    
    ggsave(filename = paste0(dir_list[f],  file_name[f], positionname(p), 
                             "intensity_change_per_cycle", '.png'),
           width = 40, height = 20, units = "cm", plot = P4)
  }
}

# RFL vs Noise analysis ---------------------------------------------------

#file path for the tsv file and create path for the tsv file
relative_path = '_fbqc_RelativeIntensityMean.tsv$'
noise_path = '_fbqc_NoiseMean.tsv$'


relative_list = c()
noise_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = relative_path, full.names = TRUE, recursive = TRUE)
  y = list.files(path =  paste0(file_path1,  e),
                 pattern = noise_path, full.names = TRUE, recursive = TRUE)
  relative_list = c(relative_list, x)
  noise_list = c(noise_list, y)
}

file_name = unlist(rm_between(relative_list, file_path1, '/', extract=TRUE))

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (i in c(1: length(file_name))){
  x = paste0(file_path2, 'SNR_plot/', file_name[i], '_SNR_plot/')
  dir_list = c(dir_list, x)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)


for (f in c(1 : length(relative_list))){
  
  dfr = read.delim(file = relative_list[f], sep = '\t', header = TRUE, 
                   stringsAsFactors = FALSE) %>% drop_na() %>% filter(Tile != "Mean")
  dfn = read.delim(file = noise_list[f], sep = '\t', header = TRUE, 
                   stringsAsFactors = FALSE) %>% drop_na() %>% filter(Tile != "Mean")
  
  for (p in unique(str_extract(dfr$Tile, "^\\w"))){
    
    dfr1 = dfr  %>% filter(str_detect(Tile, p)) %>% mutate(QC = 'RFL')
    dfn1 = dfn %>% filter(str_detect(Tile, p)) %>% mutate(QC = 'Noise')
    
    df = bind_rows(dfr1, dfn1)
    df = df %>% mutate(Cycle = as.numeric(gsub("Inc", "", Cycle))) %>%  group_by(Cycle, Base, QC) %>% 
      summarise(Value = mean(QC.Value)) %>% 
      mutate(Base=recode(Base, "A" = "G1", "T" = "G2", "G" = "R3", "C" = "R4")) %>%
      mutate(group = paste0(Base, "_", QC))
    
    datalist1 = list()
    channel = c("G1", "G2", "R3", "R4")
    for (i in c(1 :length(channel))) {
      r.vector = c()
      fitted_models = df %>% filter(Base == channel[i], Value != 0) %>% group_by(QC) %>% do(model = lm(log(Value) ~ Cycle, data = .))
      for (c in c(1: nrow(fitted_models))) r.vector[c] <- summary(fitted_models$model[[c]])$r.squared
      co_db = fitted_models %>% ungroup() %>% transmute(QC, Coef = map(model, tidy)) %>% unnest(Coef) %>% 
        filter(term == "Cycle") %>% select(QC, estimate) %>% mutate(r.squared = r.vector) %>% mutate(filter = channel[i])
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
      labs(title = paste0(file_name[f], positionname(p)," RFL vs. Noise"), x = "Cycle", y = "RFL") +
      facet_wrap(Base ~., ncol = 2, scales = "free_y")
    
    g <- ggplot_gtable(ggplot_build(P1))
    stript <- which(grepl('strip-t', g$layout$name))
    
    fills <- c( rep("#4363d8", 1), rep("#f032e6", 1), rep("#3cb44b", 1), rep("#f58231", 1)) 
    
    k <- 1
    
    for (i in stript) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    
    P3a = grid.arrange(tbl, P1, 
                       nrow=2, 
                       heights=c(2, 8))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), '_RFL_vs_noise','.png'), 
           width = 40, height = 20, units = "cm", plot = P3a)
    
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
  path2 = paste0(file_path2, 'Phasing_plot/')
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
  p_name = gsub("/phase/", "", str_extract(f, 'param_(.*)/'))
  f_name = unlist(rm_between(f, file_path1, '/', extract=TRUE))
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
  
  for (p in unique(str_extract(db$Tile, "^\\w"))){
    db1 = db %>% filter(str_detect(Tile, p))
    num.row = numrow(db1)
    
    P2 = ggplot(db1, aes(x = Cycle, y = Correlation, colour=Base, group = Base)) + 
      theme_bw() +
      geom_line(size = 0.7, alpha = 0.5) + 
      geom_point(size=0.5, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p), "Correlations Between Bases"), 
           x = "Cycle", y = "Correlation") +
      scale_x_continuous(breaks= pretty_breaks()) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text( vjust=0.6, angle = 50, size = x.font))+
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = y.font)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size= theme.font)) +
      theme(panel.spacing=unit(2, "lines")) +
      facet_wrap(Tile~., nrow = num.row)  
    
    ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), "correlation_between_bases_", 
                             'wrap_plot', '.png'), 
           width = 40, height = 20, units = "cm", plot = P2)
    
    
    P3 = ggplot(db1 %>% group_by(Cycle, Base) %>% summarise(mean_corr = mean(Correlation)), 
                aes(x = Cycle, y = mean_corr, colour=Base, group = Base)) + 
      theme_bw() +
      geom_line(size = 2, alpha = 0.5) + 
      geom_point(size=3, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p), "Correlation between bases (Mean)"),
           x = "Cycle", y = "Correlation") +
      scale_x_continuous(breaks= pretty_breaks()) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 28)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), "Correlation_between_bases", 
                             '_mean_plot', '.png'), 
           width = 40, height = 20, units = "cm", plot = P3)
    
    P3b = ggplot(db1 %>% group_by(Cycle, Base) %>% summarise(median_corr = median(Correlation)), 
                 aes(x = Cycle, y = median_corr, colour=Base, group = Base)) + 
      theme_bw() +
      geom_line(size = 2, alpha = 0.5) + 
      geom_point(size=3, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p), "Correlation between bases (Median)"), 
           x = "Cycle", y = "Correlation") +
      scale_x_continuous(breaks= pretty_breaks()) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 28)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), "Correlation_between_bases_", 
                             'median_plot', '.png'), 
           width = 40, height = 20, units = "cm", plot = P3b)
  }
}

#Cycle plot
for (f in c(1:length(cycle_list))){
  db = read.delim(file = cycle_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  for (p in unique(str_extract(db$Tile, "^\\w"))){
    db1 = db %>% filter(str_detect(Tile, p))
    
    db1$Cycle <- factor(db1$Cycle, levels=mixedsort(unique(db1$Cycle)))
    labels <- db1$Cycle[ sort(c(seq(1, nrow(db1)/length(unique(db1$Tile)), by= 20), 
                                length(unique(db1$Cycle))*4 ))]
    labels2 <- db1$Cycle[ sort(c(seq(1, nrow(db1)/length(unique(db1$Tile)), by= 40), 
                                 length(unique(db1$Cycle))*4 ))]
    
    num.row = numrow(db1)
    
    P5 = ggplot(db1 %>% filter(Cycle != "All"), aes(x = Cycle, y = Correlation, colour=Base, group = Base)) + 
      theme_bw() +
      geom_line(size = 0.7, alpha = 0.5) + 
      geom_point(size=0.5, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p), "Correlations Between Cycles"), 
           x = "Cycle", y = "Correlation") +
      scale_x_discrete(breaks=labels2) +
      scale_color_manual(values = c("A" = "#3cb44b", "T" = "#f58231", "G" = "#4363d8", "C" = "#f032e6")) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text( vjust=0.6, angle = 50, size = x.font))+
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = y.font)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size= theme.font)) +
      theme(panel.spacing=unit(2, "lines")) +
      facet_wrap(Tile~., nrow = num.row)  
    
    ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), "Correlation_between_cycles", 
                             '_wrap_plot', '.png'), 
           width = 40, height = 20, units = "cm", plot = P5)
    
    P6 = ggplot(db1 %>% group_by(Cycle, Base) %>% summarise(mean_corr = mean(Correlation)), 
                aes(x = Cycle, y = mean_corr, colour=Base, group = Base)) + 
      theme_bw() +
      geom_line(size = 2, alpha = 0.5) + 
      geom_point(size=3, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p), "Correlation between cycles (Mean)"),
           x = "Cycle", y = "Correlation") +
      scale_x_discrete(breaks=labels) +
      scale_color_manual(values = c("A" = "#3cb44b", "T" = "#f58231", "G" = "#4363d8", "C" = "#f032e6")) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20, vjust=0.6, angle = 50)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), "Correlation_between_cycles", 
                             '_mean_plot', '.png'), 
           width = 40, height = 20, units = "cm", plot = P6)
    
    
    P6b = ggplot(db1 %>% group_by(Cycle, Base) %>% summarise(median_corr = median(Correlation)), 
                 aes(x = Cycle, y = median_corr, colour=Base, group = Base)) + 
      theme_bw() +
      geom_line(size = 2, alpha = 0.5) + 
      geom_point(size=3, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p), 
                          "Correlation between cycles (Median)"), x = "Cycle", y = "Correlation") +
      scale_x_discrete(breaks=labels) +
      scale_color_manual(values = c("A" = "#3cb44b", "T" = "#f58231", "G" = "#4363d8", "C" = "#f032e6")) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20, vjust=0.6, angle = 50)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), "Correlation_between_cycles", 
                             '_median_plot', '.png'), 
           width = 40, height = 20, units = "cm", plot = P6b)
    
  }
}


#Cross_plot
for (f in c(1:length(cross2_list))){
  db = read.delim(file = cross2_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  for (p in unique(str_extract(db$Tile, "^\\w"))){
    db1 = db %>% filter(str_detect(Tile, p))
    
    db1$Cycles <- factor(db1$Cycles, levels=mixedsort(unique(db1$Cycles)))
    labels <- db1$Cycles[ sort(c(seq(1, nrow(db1)/length(unique(db1$Tile)), by= 10), length(unique(db$Cycles))*2 ))]
    labels2 <- db1$Cycles[ sort(c(seq(1, nrow(db1)/length(unique(db1$Tile)), by= 20), length(unique(db$Cycles))*2 ))]
    
    num.row = numrow(db1)
    
    P2 = ggplot(db1 %>% filter(Cycles != "All"), aes(x = Cycles, y = Percentage, colour=QC, group = QC)) +
      theme_bw() +
      geom_line(size = 0.7, alpha = 0.5) + 
      geom_point(size=0.5, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p)), 
           x = "Cycle", y = "Percentage") +
      scale_x_discrete(breaks=labels2) +
      scale_y_continuous(labels=scales::percent) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text( vjust=0.6, angle = 50, size = x.font))+
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = y.font)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size= theme.font)) +
      theme(panel.spacing=unit(2, "lines")) +
      facet_wrap(Tile~., nrow = num.row)  
    
    
    ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), '_cross2', '_wrap_plot','.png'), 
           width = 40, height = 20, units = "cm", plot = P2)
    
    
    P3 = ggplot(db1 %>% group_by(Cycles, QC) %>% summarise(mean_percent = mean(Percentage)), 
                aes(x = Cycles, y = mean_percent, colour=QC, group = QC)) + 
      theme_bw() +
      geom_line(size = 2, alpha = 0.5) + 
      geom_point(size=3, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p), "Phasing and Cross Talking (Mean)"), 
           x = "Cycles", y = "Percentage") +
      scale_x_discrete(breaks=labels) +
      scale_y_continuous(labels=scales::percent) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20, vjust=0.6, angle = 50)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p), "_cross2",  '_mean_plot', '.png'), 
           width = 40, height = 20, units = "cm", plot = P3)
    
    
    
    P3b = ggplot(db1 %>% group_by(Cycles, QC) %>% summarise(median_percent = median(Percentage)), 
                 aes(x = Cycles, y = median_percent, colour=QC, group = QC)) + 
      theme_bw() +
      geom_line(size = 2, alpha = 0.5) + 
      geom_point(size=3, alpha = 0.7) + 
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p), "Phasing and Cross Talking (Median)"), 
           x = "Cycles", y = "Percentage") +
      scale_x_discrete(breaks=labels) +
      scale_y_continuous(labels=scales::percent) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20, vjust=0.6, angle = 50)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f],   positionname(p), "_cross2",'_median_plot', '.png'), 
           width = 40, height = 20, units = "cm", plot = P3b)
  }
}


for (f in c(1:length(cross2_list))){
  db = read.delim(file = cross2_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  for (p in unique(str_extract(db$Tile, "^\\w"))){
    db1 = db %>% filter(str_detect(Tile, p)) %>% filter(Cycles != "All")
    db1$Cycles = gsub("-\\d+", "", db1$Cycles)
    db1$Cycles = as.numeric(db1$Cycles)
    
    num.row = numrow(db1)
    
    P5 = ggplot(db1 %>% group_by(Cycles, QC) %>% summarise(mean_percent = mean(Percentage)), 
                aes(x = Cycles, y = mean_percent, colour=QC, group = QC)) + 
      theme_bw() +
      geom_line(size = 1, alpha = 0.3) + 
      geom_point(size=3, alpha = 0.7) + 
      geom_smooth(method = "lm", se=FALSE, formula = y ~ x,  alpha = 0.3) +
      stat_poly_eq(formula = y ~ x, 
                   aes(label = paste(..eq.label.., sep = "~~~")), 
                   label.x = "left", label.y = "top", parse = TRUE, size = 10) +
      labs(title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p), "linear fit (mean)"), 
           x = "Cycles", y = "Percentage") +
      scale_x_continuous(breaks= pretty_breaks()) +
      scale_y_continuous(labels=scales::percent) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x=element_text(size=30), axis.text.x = element_text(size = 20)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = 28)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size=28)) +
      theme(plot.margin = margin(2, 2, 2, 2, "cm"))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], positionname(p),  "cross2_", 
                             'linear_fit_mean_plot', '.png'), 
           width = 40, height = 20, units = "cm", plot = P5)
    
    P6 = ggplot(db1, aes(x = Cycles, y = Percentage, colour=QC, group = QC)) +
      theme_bw() +
      geom_line(size = 0.7, alpha = 0.5) + 
      geom_point(size=0.5, alpha = 0.7) + 
      geom_smooth(method = "lm", se=FALSE, formula = y ~ x,  alpha = 0.3) +
      stat_poly_eq(formula = y ~ x,
                   aes(label = paste(..eq.label.., sep = "~~~")),
                   label.x.npc = "left",  parse = TRUE, size = 5) +
      labs(title = file_name[f], ' (', para_name[f], ')', positionname(p), "linear fit (mean)", 
           x = "Cycles", y = "Percentage") +
      scale_x_continuous(breaks= pretty_breaks()) +
      scale_y_continuous(labels=scales::percent) +
      theme(plot.title = element_text(hjust = 0.5, size=30)) +
      theme(axis.title.x =element_text(size=30), axis.text.x = element_text(vjust=0.6, angle = 50, 
                                                                            size = x.font))+
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=30), 
            axis.text.y = element_text(size = y.font)) +
      theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
      theme(strip.text = element_text(size= theme.font)) +
      theme(panel.spacing=unit(2, "lines")) +
      facet_wrap(Tile~., nrow= num.row)  
    
    
    ggsave(filename = paste0(dir_list[f], file_name[f],  positionname(p), "cross2_",  
                             'linear_fit_mean_plot_wrap', '.png'), 
           width = 40, height = 20, units = "cm", plot = P6)
  }
}  

# Phase percent -----------------------------------------------------------

#file path for the tsv file and create path for the tsv file
cross_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = '_phase_and_crosstalk_\\d+_cycle.tsv$', full.names = TRUE, recursive = TRUE)
  cross_list = c(cross_list, x)
}

phase.percent.list = cross_list[!grepl('_phase_and_crosstalk_2', cross_list)]
phase.percent.list2 = cross_list[grepl('_phase_and_crosstalk_2', cross_list)]
para_name = c()
for (f in phase.percent.list){
  para_name = c(para_name, gsub("/phase/", "", str_extract(f, 'param_(.*)/')))
}
file_name = unlist(rm_between(phase.percent.list, file_path1, '/', extract=TRUE))

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (i in c(1: length(file_name))){
  x = paste0(file_path2 ,'Phasing_plot/', file_name[i], '_phase_plot/', para_name[i], '/')
  dir_list = c(dir_list, x)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)

#
for (f in c(1:length(phase.percent.list))){
  df = read.delim(file = phase.percent.list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  df2 = read.delim(file = phase.percent.list2[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  for (p in unique(str_extract(df$Tile, "^\\w"))){
    dfa = df[complete.cases(df), ] %>% filter(str_detect(Tile, p)) %>% filter(Percentage != 0)
    df2a = df2[complete.cases(df2), ] %>% filter(str_detect(Tile, p)) %>% filter(Percentage != 0)
    db = bind_rows(dfa, df2a)
    
    db = db %>% filter(Cycles != "All") %>% separate(Cycles, c("Cycles", "Cycles2"))
    db$QC = gsub("Phasing\\+", "Phasing", db$QC)
    db$QC = gsub("Phasing-", "Prephasing", db$QC)
    db$Tile = gsub("mm", "", db$Tile)
    db$Cycles = as.numeric(db$Cycles)
    db1 = db %>% filter(Cycles <= 25, QC %in% c("Phasing Off-Diag", "Prephasing Off-Diag", "Phasing Prob", "Prephasing Prob")) %>% 
      group_by(QC, Tile) %>% summarise(mean_phase = mean(Percentage)) 
    db2 = db1 %>% summarise(mean_phase = mean(mean_phase))
    db2$Tile = "mean"
    db3 = bind_rows(db1, db2)
    db3$QC = gsub("Diag", "Diag (mean % of 25 cycles)", db3$QC)
    db3$QC = gsub("Prob", "Prob (mean % of 25 cycles)", db3$QC)
    
    db4 = db %>% filter(QC %in% c("Phasing Prob", "Prephasing Prob")) %>% 
      group_by(QC, Tile) %>% summarise(mean_phase = mean(Percentage)) 
    db5 = db4 %>% summarise(mean_phase = mean(mean_phase))
    db5$Tile = "mean"
    db6 = bind_rows(db4, db5)
    db6$QC = gsub("Prob", "Prob (mean % of All cycles)", db6$QC)
    
    db7 = bind_rows(db3, db6)
    
    
    db8 = db %>% group_by(QC, Cycles) %>% summarise(Percentage = mean(Percentage))
    db8$Tile = "mean"
    db9 = db %>% select(QC, Cycles, Percentage, Tile)
    db10 = bind_rows(db8, db9)
    db10$Cycles = as.numeric(db10$Cycles)
    
    P7 = ggplot(db10 %>% filter(QC %in% c("Phasing Off-Diag", "Prephasing Off-Diag", "Phasing Prob", "Prephasing Prob")),
                aes(x = Cycles, y =Percentage, col= QC, group = QC)) + 
      theme_bw() +
      geom_line(size = 1, alpha = 0.7) + 
      geom_point(size=1.5, alpha = 0.7) + 
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      labs(x = "Cycle", y = "Percentage",
           title = paste0(file_name[f], ' (', para_name[f], ')', positionname(p))) +
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
      facet_wrap(Tile ~., nrow = 2 )
    
    db7 = db7 %>% spread(Tile, mean_phase)
    db7 = db7[c(6,3,5,2, 4,1),]
    db7[, -1] = db7[, -1]
    db7 <- db7 %>% 
      mutate_if(is.numeric, round, digits = 2)
    
    f.size = tablefont(db10)
    
    tt3 <- ttheme_minimal(
      core=list(bg_params = list(blues9[1:6], col=NA),
                fg_params=list(fontface=2, fontsize = f.size)),
      colhead=list(fg_params=list(col="navyblue", fontface=4, fontsize = f.size)))
    
    tbl =  tableGrob(db7, rows =  NULL, theme = tt3)
    
    P7a = grid.arrange(tbl, P7,  
                       nrow=2, 
                       heights=c(3, 8))
    
    ggsave(filename = paste0(dir_list[f], file_name[f], '_', para_name[f], positionname(p),
                             'phase_percent','.png'), 
           width = 40, height = 20, units = "cm", plot = P7a)
    
  }
}  

# Summary report ----------------------------------------------------------

summary_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = '_summary_report', full.names = TRUE, recursive = TRUE)
  summary_list = c(summary_list, x)
}

para_name = c()
for (f in summary_list){
  para_name = c(para_name, gsub("/analyze_fastq/", "", str_extract(f, 'param_(.*)/')))
}
file_name = unlist(rm_between(summary_list, file_path1, '/', extract=TRUE))

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (i in c(1: length(file_name))){
  x = paste0(file_path2, 'Phasing_plot/', file_name[i], '_phase_plot/', para_name[i], '/')
  dir_list = c(dir_list, x)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)


for (f in c(1:length(summary_list))){
  df = read.delim(file = summary_list[f], sep = '\t', header = F, stringsAsFactors = FALSE, 
                  row.names = NULL, skip = 4)
  names(df) = df[2,]
  df = df[c(-1, -2),] %>% mutate(Category = recode(Category, "[Read Quality]" = "Read Quality", 
                                               "[Alignment]" = "Alignment"))
  colnames(df) = gsub("Category", "Run Statistics", colnames(df))
  
  a =  which(df$`Run Statistics` == "Read Quality") -1
  b = a +1
  c = which(df$`Run Statistics` == "Alignment") - 1
  d = which(df$`Run Statistics` == "Alignment")
  
  df1 = df[ c(1:a), ]
  names(df1) = gsub("^X", "",  names(df1))
  
  df2 = df[ c(b:c), ]
  colnames(df2) = gsub("Run Statistics", "Read Quality", colnames(df2))
  df2 = df2[c(-1),]
  names(df2) = gsub("^X", "",  names(df2))
  
  df3 = df[ c(d:nrow(df)), ]
  colnames(df3) = gsub("Run Statistics", "Alignment", colnames(df3))
  df3 = df3[c(-1),]
  names(df3) = gsub("^X", "",  names(df3))
  
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

# PF analysis -------------------------------------------------------------
#file path for the tsv file and create path for the tsv file
file_path2a = '_PF.tsv'

PF_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = 'Inc_?\\d?_PF.tsv$', full.names = TRUE, recursive = TRUE)
  PF_list = c(PF_list, x)
}

para_name = c()
for (f in phase.percent.list){
  para_name = c(para_name, gsub("/phase/", "", str_extract(f, 'param_(.*)/')))
}
file_name = unlist(rm_between(PF_list, file_path1, '/', extract=TRUE))

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (i in c(1: length(file_name))){
  x = paste0(file_path2 ,'PF/', file_name[i], '_PF_plot/', para_name[i], '/')
  dir_list = c(dir_list, x)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)


#Plotting 
for (f in c(1:length(PF_list))){
  
  db = read.delim(file = PF_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db$Tile2 = str_extract( db$Tile, "^\\w")
  db$Tile2 = gsub('t', 'top', db$Tile2)
  db$Tile2 = gsub('b', 'bottom', db$Tile2)
  tilelist = unique(db$Tile2)
  
  for (p in  c(1:length(tilelist))){
    
    P1 = ggplot(db %>% filter(Tile2 == tilelist[p]), aes(x = Cycle, y = Percentage, 
                                                         col=Tile, group = Tile)) + 
      theme_bw() +
      geom_line(size = 1, alpha = 0.5) + 
      geom_point(size= 2, alpha = 0.5) + 
      stat_summary(fun = mean, aes(y = Percentage, col="mean", group=1),
                   geom = 'line', size=2) +
      stat_summary(fun = mean, aes(y = Percentage, col="mean", group=1),
                   geom = 'point', size=2.5) +
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      labs(x = "Cycle", y = "%PF", 
           title = paste0(file_name[f], ' ', '(', para_name[f], ') ', 'Passing Filter'), 
           col = "Tile") +
      theme(plot.title = element_text(hjust = 0.5, size=20)) +
      theme(axis.title.x=element_text(size=20), axis.text.x = element_text(size = 20)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=20), 
            axis.text.y = element_text(size = 20)) +
      theme(legend.title = element_text(size = 20), legend.text=element_text(size=17)) +
      theme(strip.text = element_text(size=20)) +
      theme(plot.margin = margin(1, 1, 1, 1, "cm")) + 
      theme(legend.position="bottom") +
      guides(col = guide_legend(nrow = legend.row)) 
    
    ggsave(filename = paste0(dir_list[f], file_name[f], '_PF_Percent_', tilelist[p], '.png'),
           width = 40, height = 20, units = "cm", plot = P1)
    
  }
}    

# BWA analysis ------------------------------------------------------------

#file path for the tsv file and create path for the tsv file
file_path2a = '_bwa_match_tile.tsv'
file_path2b = '_perfect_match_tile.tsv'

bwa_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = file_path2a, full.names = TRUE, recursive = TRUE)
  bwa_list = c(bwa_list, x)
}

perfect_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = file_path2b, full.names = TRUE, recursive = TRUE)
  perfect_list = c(perfect_list, x)
}

para_name1 = c()
for (f in bwa_list){
  para_name1 = c(para_name1, gsub("/match/", "", str_extract(f, 'param_(.*)/')))
}
file_name1 = unlist(rm_between(bwa_list, file_path1, '/', extract=TRUE))

para_name2 = c()
for (f in perfect_list){
  para_name2 = c(para_name2, gsub("/match/", "", str_extract(f, 'param_(.*)/')))
}
file_name2 = unlist(rm_between(perfect_list,  file_path1, '/', extract=TRUE))

#creat path for plot directory, tsv directory, filename 
dir_list1 = c()
for (i in c(1: length(file_name1))){
  x = paste0( file_path2, 'BWA/', file_name1[i], '_BWA_plot/', para_name1[i], '/')
  dir_list1 = c(dir_list1, x)
}

dir_list2 = c()
for (i in c(1: length(file_name2))){
  x = paste0( file_path2, 'BWA/', file_name2[i], '_BWA_plot/', para_name2[i], '/')
  dir_list2 = c(dir_list2, x)
}

#Create directory
sapply(dir_list1, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), 
       USE.NAMES=FALSE)
sapply(dir_list2, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), 
       USE.NAMES=FALSE)


#Plotting 
for (f in c(1:length(bwa_list))){
  
  db = read.delim(file = bwa_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db$Pos = str_extract(db$Tile ,"^\\w")
  db$Pos = gsub('t', 'top', db$Pos)
  db$Pos = gsub('b', 'bottom', db$Pos)
  tilelist = unique(db$Pos)
  
  for (p in  c(1:length(tilelist))){
    
    P1 = ggplot(db %>% filter(Pos == tilelist[p]), aes(x = Cumulative.Length, y = Count, col=Tile, group = Tile)) + 
      theme_bw() +
      geom_line(size = 1, alpha = 0.7) + 
      geom_point(size=1.5, alpha = 0.7) + 
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      labs(x = "Cumulative Length", y = "Count", 
           title = paste0(file_name1[f], ' ', '(', para_name1[f], ') ', 'BWA Match'), col = "Tile") +
      theme(plot.title = element_text(hjust = 0.5, size=20)) +
      theme(axis.title.x=element_text(size=20), axis.text.x = element_text(size = 20)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=20), 
            axis.text.y = element_text(size = 20)) +
      theme(legend.title = element_text(size = 20), legend.text=element_text(size=17)) +
      theme(strip.text = element_text(size=20)) +
      theme(plot.margin = margin(1, 1, 1, 1, "cm")) + 
      theme(legend.position="bottom") +
      guides(col = guide_legend(nrow = legend.row)) 
    
    ggsave(filename = paste0(dir_list1[f], file_name1[f], '_BWA_match_' , "REF1_Total_", tilelist[p], '.png'), 
           width = 40, height = 20, units = "cm", plot = P1)
    
  }
}

#Plotting 
for (f in c(1:length(perfect_list))){
  
  db = read.delim(file = perfect_list[f], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  db$Pos = str_extract(db$Tile ,"^\\w")
  db$Pos = gsub('t', 'top', db$Pos)
  db$Pos = gsub('b', 'bottom', db$Pos)
  tilelist = unique(db$Pos)
  
  for (p in  c(1:length(tilelist))){
    
    P1 = ggplot(db %>% filter(Pos == tilelist[p]), aes(x = Cumulative.Length, y = Count, col=Tile, group = Tile)) + 
      theme_bw() +
      geom_line(size = 1, alpha = 0.7) + 
      geom_point(size=1.5, alpha = 0.7) + 
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      labs(x = "Cumulative Length", y = "Count", 
           title = paste0(file_name[f], ' ', '(', para_name[f], ') ', 'Perfect Match'), col = "Tile") +
      theme(plot.title = element_text(hjust = 0.5, size=20)) +
      theme(axis.title.x=element_text(size=20), axis.text.x = element_text(size = 20)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=20), 
            axis.text.y = element_text(size = 20)) +
      theme(legend.title = element_text(size = 20), legend.text=element_text(size=17)) +
      theme(strip.text = element_text(size=20)) +
      theme(plot.margin = margin(1, 1, 1, 1, "cm")) + 
      theme(legend.position="bottom") +
      guides(col = guide_legend(nrow = legend.row)) 
    
    ggsave(filename = paste0(dir_list2[f], file_name2[f], '_perfect_match_' , "REF1_Total_", tilelist[p], '.png'), 
           width = 40, height = 20, units = "cm", plot = P1)
    
  }
}

# Perfect decrease --------------------------------------------------------

file_path2a = '_perfect_match_tile.tsv'

perfect_list = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1,  e),
                 pattern = file_path2a, full.names = TRUE, recursive = TRUE)
  perfect_list = c(perfect_list, x)
}

#experiment_list = sapply(experiment, toString)
para_name = c()
for (f in perfect_list){
  para_name = c(para_name, gsub("/match/", "", str_extract(f, 'param_(.*)/')))
}
file_name = unlist(rm_between(perfect_list,  file_path1, '/', extract=TRUE))

#creat path for plot directory, tsv directory, filename 
dir_list = c()
for (i in c(1: length(file_name))){
  x = paste0( file_path2, 'BWA/', file_name[i], '_BWA_plot/', para_name[i], '/')
  dir_list = c(dir_list, x)
}

#Create directory
sapply(dir_list, function(x) dir.create(x, showWarnings = TRUE, recursive = TRUE, mode = "0777"), USE.NAMES=FALSE)

for (i in c(1:length(perfect_list))){
  db = read.delim(file = perfect_list[i], sep /Users/chiyuenwong/Dropbox/R_Code/Sequlite/network_E_coli_pipeline_revision3.R= '\t', header = TRUE, stringsAsFactors = FALSE)
  db1 = db %>% filter(Count != 0)
  
  for (p in unique(str_extract(db1$Tile, "^\\w"))){
    
    db1a = db1[str_detect(db1$Tile, p),] 
    
    if (max(db1a$Cumulative.Length) < 25){
      next
    }
    
    db2 = db1a %>% filter(Cumulative.Length >= 21) %>% group_by(Tile) %>% 
      mutate(percent_decrease = -(Count - lag(Count))*100/Count)
    db2 = na.omit(db2)
    
    tilelist = unique(db2$Tile)
    
    datalist2 = list()
    y = 1
    for (t in c(1: length(tilelist))){
      r.vector = c()
      db3 = db2 %>% filter(Tile == tilelist[t])
      fitted_models =  db3 %>% do(model = lm(percent_decrease ~ Cumulative.Length, data = .)) 
      co_db = fitted_models %>% ungroup() %>% transmute(Tile, Coef = map(model, tidy)) %>% unnest(Coef) %>% select(Tile, estimate)
      for (c in c(1: nrow(fitted_models))) r.vector[c] <- summary(fitted_models$model[[c]])$r.squared
      mdb = data.frame(Tile = co_db[1,1], Intercept = co_db[1,2], Slope = co_db[2,2], R.square = r.vector)
      colnames(mdb) = c("Tile", "Intercept", "Slope", "R.square")
      datalist2[[y]] = mdb
      y = y +1
    }
    
    
    db4 = bind_rows(datalist2)
    dmean = data.frame(Tile = "mean", Intercept = mean(db4$Intercept), Slope = mean(db4$Slope), R.square = mean(db4$R.square))
    db4 = bind_rows(db4, dmean)
    db4$SlopePercent = db4$Slope*100
    colnames(db4) = c("Tile", "PRD_percent",  "PRD_change_rate", "R.square", "PRD_change_rate_percent")
    db5 = db4 %>% select(Tile, PRD_percent, PRD_change_rate_percent) %>% gather(Type, Value, -Tile)
    
    P1 = ggplot(db5, aes(x = Tile, y = Value, colour=Type, group = Type)) + 
      theme_bw() +
      geom_line(size = 1.2, alpha = 0.7) + 
      geom_point(size=4, alpha = 0.7) + 
      scale_y_continuous() +
      geom_text(aes(label = round(Value, 3)),
                vjust = -0.5,
                show.legend = FALSE, size = 5) +
      labs(title = paste0(file_name[i], '_(', para_name[i], ')', positionname(p), 'perfect_read_decrease (PRD)'), 
           x = "Tile", y = "Perfect read decrease rate (%)") +
      theme(plot.title = element_text(hjust = 0.5, size=25)) +
      theme(axis.title.x=element_text(size=25), axis.text.x = element_text(size = 18, angle = 45, hjust = 0.9)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=25), 
            axis.text.y = element_text(size = 22)) +
      theme(legend.title = element_text(size = 25), legend.text = element_text(size = 20)) +
      theme(strip.text = element_text(size=25)) +
      theme(legend.position="bottom") +
      guides(col = guide_legend(nrow = 2))
    
    ggsave(filename = paste0(dir_list[i], file_name[i], '_', para_name[i], positionname(p), '_perfect_read_decrease_1', 
                             '.png'), 
           width = 40, height = 20, units = "cm", plot = P1)
    
    
    P2 = ggplot(db2, aes(x = Cumulative.Length, y = percent_decrease, colour=Tile, group = Tile)) + 
      theme_bw() +
      geom_point(size=4, alpha = 0.7) + 
      scale_y_continuous() +
      geom_smooth(method = "lm", se = FALSE) +
      labs(title = paste0(file_name[i], '_(', para_name[i], ')', positionname(p), 'perfect_match_decrease (PRD)'),
           x = "Cycles", y = "Perfect read decrease rate (%)") +
      theme(plot.title = element_text(hjust = 0.5, size=25)) +
      theme(axis.title.x=element_text(size=25), axis.text.x = element_text(size = 18)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=25), 
            axis.text.y = element_text(size = 22)) +
      theme(legend.title = element_text(size = 25), legend.text = element_text(size = 15)) +
      theme(strip.text = element_text(size=25)) +
      theme(legend.position="bottom") +
      guides(col = guide_legend(nrow = 3))
    
    db6 = db4%>% select(Tile, PRD_percent, PRD_change_rate_percent, R.square)
    rownames(db6) = db6[, 1]
    db6 = db6%>% select(PRD_change_rate_percent, PRD_percent, R.square)
    db6 = as.data.frame(t(as.matrix(db6)))
    db6 = round(db6, 3)
    db6 = setDT(db6, keep.rownames = TRUE)[]
    colnames(db6) = gsub("mm", "",  colnames(db6))
    colnames(db6) = gsub("rn", "",  colnames(db6))
    
    f.size = tablefont(db4)
    
    tt3 <- ttheme_minimal(
      core=list(bg_params = list(blues9[1:6], col=NA),
                fg_params=list(fontface=4, fontsize = f.size)),
      colhead=list(fg_params=list(col="navyblue", fontface=4, fontsize = f.size)))
    
    tbl =  tableGrob(db6, rows =  NULL, theme = tt3)
    
    P2a = grid.arrange(tbl, P2,  
                       nrow=2, 
                       heights=c(3, 8))
    
    ggsave(filename = paste0(dir_list[i], file_name[i], '_', para_name[i], positionname(p), '_perfect_read_decrease_2',
                             '.png'), width = 40, height = 20, units = "cm", plot = P2a)
  }
}  


# # Coverage analysis -------------------------------------------------------
# 
# x.font = 16
# y.font = 18
# num.col = 4
# theme.font = 16
# angle = 90
# 
# file_list = c()
# for (e in experiment){
#   path1 = paste0(file_path_coverage, e)
#   d_list = list.dirs(path1, recursive = FALSE, full.names = TRUE)
#   d_list = paste0(d_list[ grepl("param_", d_list)], "/bwa")
#   if (length(list.files(d_list)) > 0 ){
#     file_list = c(file_list, d_list)
#   }
# }
# 
# if (length(file_list) > 0){
#   
#   file_list2 = c()
#   p_name = c()
#   f_name = c()
#   for (i in file_list){
#     file_list2 = c(file_list2, list.dirs(i,  recursive = F))
#     f_name = c(f_name, list.dirs(i,  recursive = F, full.names = F))
#   }
#   
#   for (i in file_list2){
#     p_name = c(p_name, unlist(rm_between(i, str_extract(i, '_Inc_?\\d?/'), '/bwa', extract=TRUE)))
#   }
#   
#   for (f in c(1: length(file_list2))){
#     
#     filepath = list.files(path =  file_list2[f], pattern = '\\.bam$', full.names = TRUE)
#     bam.name = list.files(path =  file_list2[f], pattern = '\\.bam$', full.names = F)
#     exp = unlist(rm_between(filepath[1], file_path_coverage, '/param', extract=TRUE))
#     dir = paste0(file_path2, "Coverage/", exp, "/", p_name[f], "/", f_name[f], "/")
#     dir.create(dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
#     tilelist = c()
#     for (i in bam.name[-1]){
#       tilelist = c(tilelist, unlist(rm_between(i, str_extract(i, '_Inc_?\\d?_'), '.bam', extract=TRUE)))
#     }
#     
#     meanlist = c()
#     coverlist = c()
#     sumlist = c()
#     
#     
#     #combine plot
#     f1 = BamFile(filepath[1])
#     
#     which <- IRangesList(REF1=IRanges(1, 4583636))
#     sbp <- ScanBamParam(which = which)
#     
#     p_param <- PileupParam(distinguish_nucleotides = FALSE,
#                            distinguish_strands = FALSE,
#                            min_mapq = 0,
#                            min_nucleotide_depth= 0,
#                            min_base_quality= 0)
#     
#     res <- pileup(f1, scanBamParam = sbp, pileupParam = p_param)
#     
#     coverage.number = format(length(unique(res$pos)), big.mark=",")
#     coverage.percent = percent(length(unique(res$pos))/4583636, 0.01)
#     mean.read = format(round(sum(res$count)/4583636, 2), nsmall = 0.1)
#     IQR.value = IQR(res$count)
#     
#     tbl = tableGrob(data.frame(Coverage_breadth = coverage.number, 
#                                Coverage_breadth_percent = coverage.percent, 
#                                Coverage_depth = mean.read,
#                                Inter_quartile_range = IQR.value), rows=NULL)
#     
#     cp = ggplot(res , aes(x=count)) +
#       theme_bw() +
#       geom_histogram(bins = length(unique(res$count)), color="darkblue", fill="lightblue") +
#       labs(x = "Coverage Depth", y = "Number of Reference Bases", 
#            title = paste0(exp, " (", p_name[f], ") ", f_name[f])) +
#       theme(plot.title = element_text(hjust = 0.5, size=20)) +
#       theme(axis.title.x=element_text(size=20), axis.text.x = element_text(size = 20)) +
#       theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=20), 
#             axis.text.y = element_text(size = 15)) +
#       theme(legend.title = element_text(size = 20), legend.text=element_text(size=20)) +
#       theme(strip.text = element_text(size=20)) +
#       theme(plot.margin = margin(1, 1, 1, 1, "cm")) 
#     
#     cp2 = grid.arrange(tbl, cp, 
#                        nrow=2,
#                        heights=c(1,10))
#     
#     
#     ggsave(filename = paste0(dir, exp, '_combine_coverage', '.png'),
#            width = 40, height = 20, units = "cm", plot = cp2)
#     
#     for (i in c(2 : length(filepath))){
#       
#       f2 = BamFile(filepath[i])
#       
#       which <- IRangesList(REF1=IRanges(1, 4583636))
#       sbp <- ScanBamParam(which = which)
#       
#       p_param <- PileupParam(distinguish_nucleotides = FALSE,
#                              distinguish_strands = FALSE,
#                              min_mapq = 0,
#                              min_nucleotide_depth= 0,
#                              min_base_quality= 0)
#       
#       res <- pileup(f2, scanBamParam = sbp, pileupParam = p_param)
#       
#       sumlist = c(sumlist, sum(res$count))
#       meanlist = c(meanlist, mean(res$count))
#       coverlist = c(coverlist, length(unique(res$pos))/4583636)
#       
#     }
#     
#     df = data.frame(Tile = tilelist, Mean.coverage = meanlist, Coverage = coverlist, total.read = sumlist)
#     df$Position = str_extract(df$Tile, "^\\w")
#     
#     df$Position[df$Position == "b"] <- "bottom" 
#     df$Position[df$Position == "t"] <- "top"
#     
#     p<-ggplot(data=df %>% filter(Tile != "combine"), aes(x=Tile, y=Mean.coverage, fill=Position)) +
#       theme_bw() +
#       geom_bar(stat="identity") +
#       theme(plot.title = element_text(hjust = 0.5, size=30)) +
#       theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font, angle = angle, vjust=0.1)) +
#       theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 20), 
#             axis.text.y = element_text(hjust = 1, size = y.font)) +
#       theme(strip.text = element_text(size = theme.font)) +
#       theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
#       labs(title = paste0(exp, " (", p_name[f], ") ", f_name[f]), 
#            x = "Tile", y = "Mean Depth") +
#       theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 
#     
#     ggsave(filename = paste0(dir, exp, '_mean_depth_each_tile', '.png'),
#            width = 40, height = 20, units = "cm", plot = p)
#     
#     p1 <-ggplot(data=df %>% filter(Tile != "combine"), aes(x=Tile, y=total.read, fill=Position)) +
#       theme_bw() +
#       geom_bar(stat="identity") +
#       theme(plot.title = element_text(hjust = 0.5, size=30)) +
#       theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font, angle = angle, vjust=0.1)) +
#       theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 20), 
#             axis.text.y = element_text(hjust = 1, size = y.font)) +
#       theme(strip.text = element_text(size = theme.font)) +
#       theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
#       labs(title = paste0(exp, " (", p_name[f], ") ", f_name[f]), 
#            x = "Tile", y = "Number of Sequenced Bases \n Mapped to Reference Bases") +
#       theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 
#     
#     ggsave(filename = paste0(dir, exp, '_sequence_aligned_each_tile', '.png'),
#            width = 40, height = 20, units = "cm", plot = p1)
#     
#     
#     p2<-ggplot(data=df %>% filter(Tile != "combine"), aes(x=Tile, y=Coverage, fill=Position)) +
#       theme_bw() +
#       geom_bar(stat="identity") +
#       theme(plot.title = element_text(hjust = 0.5, size=30)) +
#       theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font, angle = angle, vjust=0.1)) +
#       theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 20), 
#             axis.text.y = element_text(hjust = 1, size = y.font)) +
#       theme(strip.text = element_text(size = theme.font)) +
#       scale_y_continuous(labels=scales::percent) +
#       theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
#       labs(title = paste0(exp, " (", p_name[f], ") ", f_name[f]), 
#            x = "Tile", y = "Coverage Breadth") +
#       theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 
#     
#     ggsave(filename = paste0(dir, exp, '_coverage_breadth_each_tile', '.png'),
#            width = 40, height = 20, units = "cm", plot = p2)
#   }
# }
# 
# 
# 
