# This script is for assigning ChIP peaks to nearest genes (within 100kb of TSS of nearest genes)
library(dplyr)

# set path for Mac
# dir <- "Dropbox (OIST)/Ishikawa Unit/Tsunghan/180604_Th1_ChIP"

# set path for Windows
 dir <- "E:/Dropbox (OIST)/Ishikawa Unit/Tsunghan/180604_Th1_ChIP/"

# dir <- "/Users/Tsunghan/Dropbox (OIST)/Ishikawa Unit/Tsunghan/180604_Th1_ChIP/annotated"

# functions
ChIP_filter_fun <- function(x) {
  x <- x %>%
    filter(abs(as.numeric(as.character(Distance.to.TSS))) < 100000) %>%
    select(Gene.Name) %>%
    distinct(Gene.Name)
}

JunB <- read.delim(file.path(dir,"annotated",paste("Th1_Th17_Treg_JunB",".bed",sep="")),
                   sep="\t",header=TRUE,fill=TRUE,dec=".")

BATF <- read.delim(file.path(dir,"annotated",paste("Th1_Th17_Treg_BATF",".bed",sep="")),
                   sep="\t",header=TRUE,fill=TRUE,dec=".")

IRF4 <- read.delim(file.path(dir,"annotated",paste("Th1_Th17_Treg_IRF4",".bed",sep="")),
                   sep="\t",header=TRUE,fill=TRUE,dec=".")

overlapped <- read.delim(file.path(dir,"annotated",paste("111_Th1_Th17_Treg_JunB_BATF_IRF4_filtered_peaks",".bed",sep="")),
                   sep="\t",header=TRUE,fill=TRUE,dec=".")


ChIP_list <- list(JunB, BATF, IRF4, overlapped)
ChIP_filtered <- lapply(ChIP_list,ChIP_filter_fun)

write.csv(ChIP_filtered[[1]],file=file.path(dir,paste("Th1_Th17_Treg_JunB_gene_list",".csv",sep="")))
write.csv(ChIP_filtered[[2]],file=file.path(dir,paste("Th1_Th17_Treg_BATF_gene_list",".csv",sep="")))
write.csv(ChIP_filtered[[3]],file=file.path(dir,paste("Th1_Th17_Treg_IRF4_gene_list",".csv",sep="")))
write.csv(ChIP_filtered[[4]],file=file.path(dir,paste("Th1_Th17_Treg_JunB_BATF_IRF4_gene_list",".csv",sep="")))
