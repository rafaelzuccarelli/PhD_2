library(car)
library(ggplot2)

#Store arguments from terminal
args <- commandArgs(trailingOnly=TRUE)

#data <- read.csv("/home/rafael/Bioinfo/R/Productivity/Plant/ALL/Productivity_all.csv", header = TRUE, comment.char="#")
#get the file name and path from terminal 
wd<-getwd()
slash<-'/'
#the paste function join together the working directory (wd) a slash(variable created above) and the csv file name given in terminal
path_file <- paste(wd,slash,args[1],sep= "",collapse=NULL)
#after correctly store the full path name created above stores the file content in data
data <- read.csv(path_file, header = TRUE, comment.char="#")

#slicing the genotypes using the name of the genotypes in the the first collumn
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

#column of the file that are going to be analysed. Creates a file for output with the same name plus .csv. Starts with 2 because 1 is the samples names.
#This is the part that recursevely parse over all collumns.

n1 <- 2 #firt collumn with data
colname <-colnames(data[n1])
file_name <- paste(colname, ".csv",sep= "",collapse=NULL)
output <- file(file_name, "w")
numcol <- ncol(data)#number of collumns in the file.

while(n1<=numcol){
  colname <-colnames(data[n1])#name of the collumn (header)
  file_name <- paste(colname, ".csv",sep= "",collapse=NULL)#file name created using the header of collumn.
  output <- file(file_name, "w")
  #slicing the measurements, one collumn to each genotype (already sliced above), n1 correspond to the collumn number that are going to increment in the while loop.
  #This loop will repeat this until the number of the collumns in the file ends.
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
  #sample names
  NG01 <- WT[1,1]
  NG02 <- RNAi_35S_L2[1,1]
  NG03 <- RNAi_35S_L5[1,1]
  NG04 <- OE_35S_L1[1,1] 
  NG05 <- OE_35S_L13[1,1]
  NG06 <- OE_35S_L37[1,1] 
  NG07 <- RNAi_PPC2_L5A[1,1]
  NG08 <- RNAi_PPC2_L5B[1,1]
  NG09 <- OE_PPC2_L2B[1,1]
  NG10 <- OE_PPC2_L20B[1,1]
  NG11 <- OE_PPC2_L21[1,1] 
  NG12 <- OE_PG_L1[1,1]
  NG13 <- OE_PG_L10[1,1]
  
  #slice list whithout WT (because all samples will be compared with WT)
  SL <- list(RG02,RG03,RG04,RG05,RG06,RG07,RG08,RG09,RG10,RG11,RG12,RG13)
  NL <- list(NG02,NG03,NG04,NG05,NG06,NG07,NG08,NG09,NG10,NG11,NG12,NG13)
  CL <- list(1,2,3,4,5,6,7,8,9,10,11,12)
  
  txt1a<-'WT'
  txt2a<-toString(mean(RG01, na.rm=TRUE))
  txt3a<- ','
  semRG01<-sd(RG01,na.rm=TRUE)/sqrt(length(RG01))
  txt4a<-toString(semRG01)
  txt5a<-'-'
  
  txt6a<-paste(txt1a,txt3a,txt2a,txt3a,txt4a,txt3a,txt5a,txt5a,txt5a, sep= "", collapse=NULL)
  
  txt7a<-'Variable'
  txt8a<-'Mean'
  txt9a<-'SE'
  txt10a<-'p_value'
  txt11a<-'test'
  txt12a<-toString(colname)
  txt13a<-paste(txt7a,txt3a,txt8a,txt3a,txt9a,txt3a,txt10a,txt3a,txt11a,txt3a,txt12a,sep= "", collapse=NULL)
  
  write(txt13a, output,append = TRUE)
  write(txt6a, output,append = TRUE)
  
  
  #shapiro–Wilk normality test, if p > 0.05 the sample is normal 
  for(i in CL){
    shapiro.test(RG01) -> WT_S
    shapiro.test(unlist(SL[i])) -> RGX_S
    #constructing a consistent .csv file showing mean, SE, p-value and test (T or W)
    if(WT_S$p.value > 0.05 & RGX_S$p.value > 0.05){ #test if both are normal, if yes do the t.test, if not do the wilcox test)
      #all the bellow code is to form a .csv file easy to export to origins program
      t.test(RG01,unlist(SL[i]),conf.level = 0.95, var.equal = FALSE,na.rm=TRUE) -> t_test_result
      #print(t_test_result)
      txt1b<-unlist(NL[i])
      txt2b<-toString(mean(unlist(SL[i]), na.rm=TRUE))
      txt3b<- ','
      semi<-sd(unlist(SL[i]), na.rm=TRUE)/sqrt(length(unlist(SL[i])))
      txt4b<-toString(semi)
      txt5b<-toString(t_test_result$p.value)
      txt6b<-'t_test'
      txt7b<-paste(txt1b,txt3b,txt2b,txt3b,txt4b,txt3b,txt5b,txt3b,txt6b, sep= "",collapse=NULL)
      write(txt7b,output,append = TRUE)
    }else{
      wilcox.test(RG01,unlist(SL[i])) -> wilcox_results
      #print(wilcox_results)
      txt1b<-unlist(NL[i])
      txt2b<-toString(mean(unlist(SL[i]), na.rm=TRUE))
      txt3b<- ','
      semi<-sd(unlist(SL[i]), na.rm=TRUE)/sqrt(length(unlist(SL[i])))
      txt4b<-toString(semi)
      txt5b<-toString(wilcox_results$p.value)
      txt6b<-'Wilcox'
      txt7b<-paste(txt1b,txt3b,txt2b,txt3b,txt4b,txt3b,txt5b,txt3b,txt6b, sep= "",collapse=NULL)
      write(txt7b, output, append = TRUE)
    }  
  }
  #boxplot
  data_frame <- data.frame(data[1],data[n1])
  x<-colnames(data[1])
  p <- ggplot(data_frame, aes(x= !!ensym(x), y= !!ensym(colname), fill= !!ensym(x))) + geom_boxplot(show.legend=FALSE)+
    scale_fill_manual(values=c("#999999","#999999","#999999","#999999","#999999","#999999","#999999","#999999","#999999","#999999","#999999","#999999","#666666"))+
    scale_x_discrete(limits=c("WT","35S_RNAi_L2","35S_RNAi_L5","35S_OE_L1","35S_OE_L13","35S_OE_L37","PPC2_RNAi_L5A","PPC2_RNAi_L5B","PPC2_OE_L2B","PPC2_OE_L20B","PPC2_OE_L21","PG_OE_L1","PG_OE_L10")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.background = element_rect(fill = "white", colour = NA))
  boxplot_name<-paste(colname, ".png",sep= "",collapse=NULL)
  ggsave(boxplot_name, width = 20, height = 10, units = "cm")
  
  n1 <- n1+1 #add one to the collumn numbe, therefore jump to next collumn
}

