rm(list=ls(all=TRUE))

library(survival)
library(MASS)
library(parallel)

#1940-2018
lalo <- read.table("latlonERA5.txt",head=FALSE)
c=1
cc=1
# water year 10.1-9.30


data.set<-read.table("starwin1940ERA5.txt",head=FALSE)#366
colnames(data.set) <- c("STAR")
STOP<-read.table("stopwin1940ERA5.txt",head=FALSE)

namens<-paste("pdo1940-2018_ndjfm_smoo7daywinter.txt",sep="")
PDO <- read.table (namens, header=F, sep=",")
PDO=scale(PDO)
data.set$PDO  <- PDO


#Reading the data

coeff.coxPDO=matrix(-9999,lalo[c,1],lalo[c,2])
pval.coxPDO=matrix(-9999,lalo[c,1],lalo[c,2])


 
for(i in 120:(lalo[cc,1]+120-1))
{
  
  na1<-paste("/Users/zhiqiyang/Documents/A-ucsd-postdoc/ERA5winteronly/arevent194011_20183_lat",i,"indp_nov2mar_alldomain.txt",sep="")
  evall <- read.table (na1, header=F, sep="")

  
  for(j in 1:(lalo[cc,2]))
  {
    evet <- evall[,j]
    evet<-evet[-c(1:3)]
    
    
    data.set$STOP <- STOP
    data.set$EVENT <- evet
    
    
    if(sum(data.set[,4])>=100)
    {
      #Model selection
      mod.b=coxph(Surv(STAR,STOP$V1,EVENT) ~ 1,data=data.set)
      mod.m1=coxph(Surv(STAR,STOP$V1,EVENT) ~ PDO,data=data.set)
      
      
      
      
      e1=extractAIC(mod.m1)[2]
      
      eb=extractAIC(mod.b)[2]
      
      
      
      if(e1<eb)
      {
        
        
          fmod=mod.m1
          coeff.coxPDO[i-120+1,j]=coef(summary(fmod))[,1]
          pval.coxPDO[i-120+1,j]=coef(summary(fmod))[,5]
        
  
        
        
      }
    }
  }
}






namePDOc <-paste("ARregion_coeffPDO1940-2018indp_nov2mar.txt",sep="")  
namePDOp <-paste("ARregion_pvalPDO1940-2018indp_nov2mar.txt",sep="") 
write.table(coeff.coxPDO, file = namePDOc, col.names = F,row.names = F, quote = F)
write.table(pval.coxPDO, file = namePDOp, col.names = F,row.names = F, quote = F)




