# This script is used for illustrating a scatter plot from de novo motif analysis by Homer 

library(dplyr)
library(ggrepel)
library(ggplot2)
library(scales)
library(tidyr)

# set dir
# dir = "Dropbox (OIST)/Ishikawa Unit/Tsunghan/180604_Th1_ChIP/filtered_motif/"
dir <- "E:/Dropbox (OIST)/Ishikawa Unit/Tsunghan/manuscript_JI/BETA_Th1_Th17_Treg/Raw"
# dir <- "/Users/Tsunghan/Dropbox (OIST)/Ishikawa Unit/Tsunghan/180604_Th1_ChIP/filtered_motif"

data_modified_fun <- function(x) {
  x <- x %>%
    select(Motif.Name,Log.P.value,X..of.Target.Sequences.with.Motif) %>%
    separate(Motif.Name,c("Motif","TF")) %>%
    mutate(percentages = gsub('%', '', x$X..of.Target.Sequences.with.Motif))
}


# read data and basic wranggling
AP_1 <- read.delim(file.path(dir,"111_Th1_Th17_Treg_JunB_BATF_IRF4_filtered_peaks.bed","knownResults.txt"),header=TRUE,fill=TRUE,sep='\t')

motif_list<-list(AP_1)
motif_list_modified <- lapply(motif_list,data_modified_fun)

data <- motif_list_modified[[1]]
data$motiflabels <- ifelse(-(data$Log.P.value) > 200, TRUE, FALSE)
ggplot(data, aes(x=as.double(percentages), y=-(Log.P.value))) +
  xlim(0,100) +
  ylim(0,4500) +
  theme_classic() +
  geom_point() +
  geom_text_repel(label = ifelse(data$motiflabels == TRUE, as.character(data$Motif),""), size = 3, box.padding = 0.1)
