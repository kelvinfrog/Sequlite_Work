library(wesanderson)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(stringr)
library(gtools)
library(ggpmisc)
library(broom)
library(uniftest)
library(ggpubr)
library(gridExtra)
library(Rsamtools)
library(qdapRegex)


x.font = 16
y.font = 18
num.col = 4
theme.font = 16
angle = 90

experiment  = c("67x_Inc")
file_path1 = 'C:/Users/kelvin/Documents/Bioinformatics_analysis/'
exp = c()
for (e in experiment) {
  x = list.files(path =  paste0(file_path1), pattern = e, full.names = TRUE, recursive = F)
  exp = c(exp, x)
}

experiment = experiment[str_detect(experiment, "Inc")]


#file path for the tsv file and create path for the tsv file
file_path1 = 'C:/Users/kelvin/Documents/Bioinformatics_analysis/'

file_list = c()
for (e in experiment){
  path1 = paste0(file_path1, "CRT", e)
  d_list = list.dirs(path1, recursive = FALSE, full.names = TRUE)
  d_list = d_list[ grepl("param_", d_list)]
  file_list = c(file_list, d_list)
}

file_list = unlist(lapply(file_list, function(x) paste0(x, "/bwa")))

file_list2 = c()
p_name = c()
f_name = c()
for (i in file_list){
  file_list2 = c(file_list2, list.dirs(i,  recursive = F))
  f_name = c(f_name, list.dirs(i,  recursive = F, full.names = F))
}

for (i in file_list2){
  p_name = c(p_name, unlist(rm_between(i, '_Inc/', '/bwa', extract=TRUE)))
}


for (f in c(1: length(file_list2))){
  

  filepath = list.files(path =  file_list2[f], pattern = '\\.bam$', full.names = TRUE)
  
  exp = unlist(rm_between(filepath[1], '_analysis/', '/param', extract=TRUE))
 
  dir = paste0("C:/Users/kelvin/Documents/Figures/Coverage/", exp, "/", p_name[f], "/", f_name[f], "/")
  
  dir.create(dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  
  tilelist = unlist(lapply(filepath, str_extract, "\\w?\\d\\d\\.?\\d?\\d?mm"))[-1]
  
  meanlist = c()
  coverlist = c()
  sumlist = c()
  
  #combine plot
  f1 = BamFile(filepath[1])
  
  which <- IRangesList(REF1=IRanges(1, 4583636))
  sbp <- ScanBamParam(which = which)
  
  p_param <- PileupParam(distinguish_nucleotides = FALSE,
                         distinguish_strands = FALSE,
                         min_mapq = 0,
                         min_nucleotide_depth= 0,
                         min_base_quality= 0)
  
  res <- pileup(f1, scanBamParam = sbp, pileupParam = p_param)
  
  coverage.number = format(length(unique(res$pos)), big.mark=",")
  coverage.percent = percent(length(unique(res$pos))/4583636, 0.01)
  mean.read = format(round(sum(res$count)/4583636, 2), nsmall = 0.1)
  IQR.value = IQR(res$count)
  
  tbl = tableGrob(data.frame(Coverage_breadth = coverage.number, 
                             Coverage_breadth_percent = coverage.percent, 
                             Coverage_depth = mean.read,
                             Inter_quartile_range = IQR.value), rows=NULL)
  
  cp = ggplot(res , aes(x=count)) +
    theme_bw() +
    geom_histogram(bins = length(unique(res$count)), color="darkblue", fill="lightblue") +
    labs(x = "Coverage Depth", y = "Number of Reference Bases", 
         title = paste0(exp, " (", p_name[f], ") ", f_name[f])) +
    theme(plot.title = element_text(hjust = 0.5, size=20)) +
    theme(axis.title.x=element_text(size=20), axis.text.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=20), 
          axis.text.y = element_text(size = 15)) +
    theme(legend.title = element_text(size = 20), legend.text=element_text(size=20)) +
    theme(strip.text = element_text(size=20)) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) 
  
  cp2 = grid.arrange(tbl, cp, 
               nrow=2,
               heights=c(1,10))
  
  
  ggsave(filename = paste0(dir, exp, '_combine_coverage', '.png'),
         width = 40, height = 20, units = "cm", plot = cp2)

  
  for (i in c(2 : length(filepath))){
    
    f2 = BamFile(filepath[i])
    
    which <- IRangesList(REF1=IRanges(1, 4583636))
    sbp <- ScanBamParam(which = which)
    
    p_param <- PileupParam(distinguish_nucleotides = FALSE,
                           distinguish_strands = FALSE,
                           min_mapq = 0,
                           min_nucleotide_depth= 0,
                           min_base_quality= 0)
    
    res <- pileup(f2, scanBamParam = sbp, pileupParam = p_param)
    
    sumlist = c(sumlist, sum(res$count))
    meanlist = c(meanlist, mean(res$count))
    coverlist = c(coverlist, length(unique(res$pos))/4583636)
    
  }
  
  
  
  df = data.frame(Tile = tilelist, Mean.coverage = meanlist, Coverage = coverlist, total.read = sumlist)
  df$Position = str_extract(df$Tile, "^\\w")
  
  
  df$Position[df$Position == "b"] <- "bottom" 
  df$Position[df$Position == "t"] <- "top"
  
  p<-ggplot(data=df %>% filter(Tile != "combine"), aes(x=Tile, y=Mean.coverage, fill=Position)) +
    theme_bw() +
    geom_bar(stat="identity") +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font, angle = angle, vjust=0.1)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 20), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size = theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(exp, " (", p_name[f], ") ", f_name[f]), 
         x = "Tile", y = "Mean Depth") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 
  
  ggsave(filename = paste0(dir, exp, '_mean_depth_each_tile', '.png'),
         width = 40, height = 20, units = "cm", plot = p)
  #print(p)
  
  p1 <-ggplot(data=df %>% filter(Tile != "combine"), aes(x=Tile, y=total.read, fill=Position)) +
    theme_bw() +
    geom_bar(stat="identity") +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font, angle = angle, vjust=0.1)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 20), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size = theme.font)) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(exp, " (", p_name[f], ") ", f_name[f]), 
         x = "Tile", y = "Number of Sequenced Bases \n Mapped to Reference Bases") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 
  #print(p1)
  
  ggsave(filename = paste0(dir, exp, '_sequence_aligned_each_tile', '.png'),
         width = 40, height = 20, units = "cm", plot = p1)
  
  
  p2<-ggplot(data=df %>% filter(Tile != "combine"), aes(x=Tile, y=Coverage, fill=Position)) +
    theme_bw() +
    geom_bar(stat="identity") +
    theme(plot.title = element_text(hjust = 0.5, size=30)) +
    theme(axis.title.x=element_text(size=28), axis.text.x = element_text(size = x.font, angle = angle, vjust=0.1)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 20), 
          axis.text.y = element_text(hjust = 1, size = y.font)) +
    theme(strip.text = element_text(size = theme.font)) +
    scale_y_continuous(labels=scales::percent) +
    theme(legend.title = element_text(size = 30), legend.text=element_text(size=28)) +
    labs(title = paste0(exp, " (", p_name[f], ") ", f_name[f]), 
         x = "Tile", y = "Coverage Breadth") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) 
  #print(p2)
  
  ggsave(filename = paste0(dir, exp, '_coverage_breadth_each_tile', '.png'),
         width = 40, height = 20, units = "cm", plot = p2)
  
}























