library(wesanderson)
library(ggplot2)
library(plyr)
library(tidyr)
library(dplyr)

#list the directory contains the tsv file

path = 'C:/Users/kelvin/Bioinformatics/321_pixel_raw_file/Pixel_Q123_union_file/'
save_path = 'C:/Users/kelvin/Bioinformatics/321_pixel_raw_file/Pixel_Q123_plot/'
file_description = '_PixelQ123'


pixel_raw_table <- function(folder_path){
  
  #make file path list  
  path_list <- function(folder_path){
    
    filezz = list.files(path)
    pathz = c()
    
    for (i in filezz)
    {
      pathi = paste0(path, "/", i)
      pathz = append(pathi, pathz)
    }  
    return(pathz)
  }
  pathlist = path_list(path)
  
  
  #make graph heading list
  title_head <- function(folder_path){
    
    LL = c()
    filez <- list.files(folder_path)
    
    for (i in filez)
    {
      a = sub('_.*', '', i) 
      LL <- append(a, LL)
    }
    return(LL)
  }
  headlist = title_head(path)
  
  #start the for loop to make graphs
  for (i in c(1:length(pathlist)))
  {
    
  #read the table  
  db = read.table(file = pathlist[i], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  
  #reorganize the table
  db = db[, -2]
  db = db %>% gather(Filter, FL, -Cycle, -QC)
  db$Cycle <- gsub('Inc', '', db$Cycle)
  db$Cycle <- as.integer(as.character(db$Cycle))
  
  P = ggplot(db, aes(x = Cycle, y = FL, group = Filter, colour=Filter)) + 
    theme(panel.background = element_rect(fill = 'gray95')) +
    geom_line(size = 1.2, alpha = 0.7) + 
    geom_point(size=4, alpha = 0.7) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = length(unique(db$Cycle)))) +
    scale_color_manual(values = c("G1" = "#3cb44b", "G2" = "#f58231", "R3" = "#4363d8", "R4" = "#f032e6")) +
    theme(axis.title.x=element_text(size=20)) +
    theme(axis.text.x = element_text( size = 15)) +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 20)) +
    theme(axis.text.y = element_text(hjust = 1, size = 15)) +
    #scale_color_brewer(palette="Accent") +
    facet_grid(. ~ QC) +
    ggtitle(paste(headlist[i],"-", "Pixel Q123")) +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  print(P)
  #ggsave(filename = paste0(save_path,headlist[i],file_description,".png"), plot = P)
  
  }
}

pixel_raw_table(path)

  
  