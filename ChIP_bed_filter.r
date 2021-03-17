# This script is for filtering ChIPseq peaks for further analysis

library(dplyr)

# set path for Mac
# dir <- "Dropbox (OIST)/Ishikawa Unit/Tsunghan/180604_Th1_ChIP/"

# set path for Windows
dir <- "/Users/Tsunghan//Dropbox (OIST)/Ishikawa Unit/Tsunghan/180604_Th1_ChIP/"
# dir <- "E:/Dropbox (OIST)/Ishikawa Unit/Tsunghan/180604_Th1_ChIP/"

# This function is for selecting peaks located within 100kb of TSS of nearest gene bodies and stored in format for Homer.
ChIP_filter_fun <- function(x) {
  x <- x %>%
    filter(abs(as.numeric(as.character(Distance.to.TSS))) < 100000)
  x <- x[,c(2,3,4,1,7,5)]
}

# peak.score == normalized tag counts

JunB <- read.delim(file.path(dir,"annotated",paste("Th1_JunB_WT_1",".txt",sep="")),
                   sep="\t",header=TRUE,fill=TRUE,dec=".")

BATF <- read.delim(file.path(dir,"annotated",paste("Th1_BATF_WT_1",".txt",sep="")),
                   sep="\t",header=TRUE,fill=TRUE,dec=".")

IRF4 <- read.delim(file.path(dir,"annotated",paste("Th1_IRF4_WT_1",".txt",sep="")),
                   sep="\t",header=TRUE,fill=TRUE,dec=".")

BATF_KO <- read.delim(file.path(dir,"annotated",paste("Th1_BATF_KO_1",".txt",sep="")),
                      sep="\t",header=TRUE,fill=TRUE,dec=".")


IRF4_KO <- read.delim(file.path(dir,"annotated",paste("Th1_IRF4_KO_1",".txt",sep="")),
                      sep="\t",header=TRUE,fill=TRUE,dec=".")

ChIP_list <- list(JunB, BATF, IRF4, BATF_KO, IRF4_KO)
ChIP_filtered <- lapply(ChIP_list,ChIP_filter_fun)

write.table(ChIP_filtered[[1]],file=file.path(dir,"filtered_peaks",paste("Th1_JunB_filtered_peaks",".bed",sep="")),sep="\t",row.names=F,col.names=F,quote=F)
write.table(ChIP_filtered[[2]],file=file.path(dir,"filtered_peaks",paste("Th1_BATF_filtered_peaks",".bed",sep="")),sep="\t",row.names=F,col.names=F,quote=F)
write.table(ChIP_filtered[[3]],file=file.path(dir,"filtered_peaks",paste("Th1_IRF4_filtered_peaks",".bed",sep="")),sep="\t",row.names=F,col.names=F,quote=F)
write.table(ChIP_filtered[[4]],file=file.path(dir,"filtered_peaks",paste("Th1_BATF_KO_filtered_peaks",".bed",sep="")),sep="\t",row.names=F,col.names=F,quote=F)
write.table(ChIP_filtered[[5]],file=file.path(dir,"filtered_peaks",paste("Th1_IRF4_KO_filtered_peaks",".bed",sep="")),sep="\t",row.names=F,col.names=F,quote=F)




