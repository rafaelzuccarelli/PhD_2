library(car)

#Store arguments from terminal
args <- commandArgs(trailingOnly=TRUE)
#data <- read.csv("/home/rafael/Bioinfo/R/Productivity/Fruit/Brix/plant.csv", header = TRUE, comment.char="#")
#get the file name and path from terminal
wd<-getwd()
slash<-'/'
path_file <- paste(wd,slash,args[1],sep= "",collapse=NULL)
data <- read.csv(path_file, header = TRUE, comment.char="#")

#slicing the genotypes
WT <- subset(data, Genotype == "WT")
RNAi_35S_L2 <- subset(data, Genotype == "35S_RNAi_L2")
RNAi_35S_L5 <- subset(data, Genotype == "35S_RNAi_L5")
OE_35S_L1 <- subset(data, Genotype == "35S_OE_L1")
OE_35S_L13 <- subset(data, Genotype == "35S_OE_L13")
OE_35S_L37 <- subset(data, Genotype == "35S_OE_L37")
RNAi_PPC2_L5A <- subset(data, Genotype == "PPC2_RNAi_L5A")
RNAi_PPC2_L5B <- subset(data, Genotype == "PPC2_RNAi_L5B")
OE_PPC2_L2B <- subset(data, Genotype == "PPC2_OE_L2B")
OE_PPC2_L20B <- subset(data, Genotype == "PPC2_OE_L20B")
OE_PPC2_L21 <- subset(data, Genotype == "PPC2_OE_L21")
OE_PG_L1 <- subset(data, Genotype == "PG_OE_L1")
OE_PG_L10 <- subset(data, Genotype == "PG_OE_L10")

#column of the file that are going to be analysed. Creates a file for output with the same name plus .csv

n1 <- 2
colname <-colnames(data[n1])
file_name <- paste(colname, ".csv",sep= "",collapse=NULL)
output <- file(file_name, "w")

numcol <- ncol(data)
while(n1<=numcol){
  colname <-colnames(data[n1])
  file_name <- paste(colname, ".csv",sep= "",collapse=NULL)
  output <- file(file_name, "w")
  #slicing the measurements 
  RG01 <- WT[,n1]
  RG02 <- RNAi_35S_L2[,n1]
  RG03 <- RNAi_35S_L5[,n1]
  RG04 <- OE_35S_L1[,n1] 
  RG05 <- OE_35S_L13[,n1]
  RG06 <- OE_35S_L37[,n1] 
  RG07 <- RNAi_PPC2_L5A[,n1]
  RG08 <- RNAi_PPC2_L5B[,n1]
  RG09 <- OE_PPC2_L2B[,n1]
  RG10 <- OE_PPC2_L20B[,n1]
  RG11 <- OE_PPC2_L21[,n1] 
  RG12 <- OE_PG_L1[,n1]
  RG13 <- OE_PG_L10[,n1]
  
  #slice list whithout WT
  SL <- list(RG02,RG03,RG04,RG05,RG06,RG07,RG08,RG09,RG10,RG11,RG12,RG13)
  
  #shapiro–Wilk normality test, if p > 0.05 the sample is normal 
  for(i in SL){
    shapiro.test(RG01) -> WT_S
    shapiro.test(i) -> RGX_S
    if(WT_S$p.value > 0.05 & RGX_S$p.value > 0.05){ #test if both are normal, if yes do the t.test, if not do the wilcox test)
      t.test(RG01,i,conf.level = 0.95, var.equal = FALSE) -> t_test_result
      #print(t_test_result)
      #constructing a consistent .csv file showing mean, SE, p-value and test (T or W)
      txt1a<-'Mean_SE WT = '
      txt2a<-toString(mean(RG01, na.rm=TRUE))
      txt3a<- ','
      semRG01<-sd(RG01,na.rm=TRUE)/sqrt(length(RG01))
      txt4a<-toString(semRG01)
      txt5a<-'T'
      txt6a<-paste(txt1a, txt3a, txt2a, txt3a, txt4a, txt3a, txt5a, sep= "",collapse=NULL)
      write(txt6a, output,append = TRUE)
      txt1b<-'Mean_SE i = '
      txt2b<-toString(mean(i, na.rm=TRUE))
      txt3b<- ','
      semi<-sd(i, na.rm=TRUE)/sqrt(length(i))
      txt4b<-toString(semi)
      txt5b<-toString(t_test_result$p.value)
      txt6b<-'T'
      txt7b<-paste(txt1b, txt3b, txt2b, txt3b, txt4b, txt3b, txt5b, txt3b, txt6b, sep= "",collapse=NULL)
      write(txt7b,output,append = TRUE)
    }else{
      wilcox.test(RG01,i) -> wilcox_results
      #print(wilcox_results)
      txt1a<-'Mean_SE WT = '
      txt2a<-toString(mean(RG01, na.rm=TRUE))
      txt3a<- ','
      semRG01<-sd(RG01, na.rm=TRUE)/sqrt(length(RG01))
      txt4a<-toString(semRG01)
      txt5a<-'W'
      txt6a<-paste(txt1a, txt3a, txt2a, txt3a, txt4a, txt3a, txt5a, sep= "",collapse=NULL)
      write(txt6a, output, append = TRUE)
      txt1b<-'Mean_SE i = '
      txt2b<-toString(mean(i, na.rm=TRUE))
      txt3b<- ','
      semi<-sd(i, na.rm=TRUE)/sqrt(length(i))
      txt4b<-toString(semi)
      txt5b<-toString(wilcox_results$p.value)
      txt6b<-'W'
      txt7b<-paste(txt1b, txt3b, txt2b, txt3b, txt4b, txt3b, txt5b, txt3b, txt6b, sep= "",collapse=NULL)
      write(txt7b, output, append = TRUE)
    }  
  } 
  n1 <- n1+1
}

