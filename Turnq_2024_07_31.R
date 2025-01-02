#Written by James Fordyce 06.03.2024

library(hillR);library(vegan);library(reshape2);library(ggplot2);library(vegetarian)

#functions
turnq<-function(mat1,mat2,name.mat1="mat1",name.mat2="mat2",cormethod="pearson",qgrain=1,seed=2){
  
  qz1<-round(seq(from=0, to=5, by=qgrain),digits=1)
  qz2<-round(seq(from=0, to=5, by=qgrain),digits=1)
  
  resm<-matrix(ncol=length(qz1),nrow=length(qz2)) #mantel
  colnames(resm)<-paste(name.mat1,qz2,sep="_")
  rownames(resm)<-paste(name.mat2,qz1,sep="_")
  
  resp<-matrix(ncol=length(qz1),nrow=length(qz2)) #procrustes
  colnames(resp)<-paste(name.mat1,qz2,sep="_")
  rownames(resp)<-paste(name.mat2,qz1,sep="_")
  
  resm.p<-matrix(ncol=length(qz1),nrow=length(qz2))
  colnames(resm.p)<-paste(name.mat1,qz2,sep="_")
  rownames(resm.p)<-paste(name.mat2,qz1,sep="_")
  
  resp.p<-matrix(ncol=length(qz1),nrow=length(qz2))
  colnames(resp.p)<-paste(name.mat1,qz2,sep="_")
  rownames(resp.p)<-paste(name.mat2,qz1,sep="_")
  
  
  set.seed(seed)
  
  for(i in 1:length(qz1)){
    for(j in 1:length(qz2)){
      #hill_taxa_parti_pairwise(mat1,q=qz1[i],output="matrix",pairs="full")$TD_beta-1
      
      a<-as.dist(hill_taxa_parti_pairwise(mat1,q=qz1[i],output="matrix",pairs="full")$TD_beta-1,diag=TRUE)
      b<-as.dist(hill_taxa_parti_pairwise(mat2,q=qz2[j],output="matrix",pairs="full")$TD_beta-1,diag=TRUE)
      
      corr<-mantel(a,b,method=cormethod)$statistic
      corr.p<-mantel(a,b,method=cormethod)$signif
      
      procor<-protest(scores(dbrda(a~1,add=TRUE)),scores(dbrda(b~1,add=TRUE)))[[6]]
      procor.p<-protest(scores(dbrda(a~1,add=TRUE)),scores(dbrda(b~1,add=TRUE)))[[13]]

      
      resm[j,i]<-corr
      resp[j,i]<-procor
      resm.p[j,i]<-corr.p
      resp.p[j,i]<-procor.p
    }
  }
  res<-list(mantel.cor=resm,mantel.p=resm.p,procrusts.cor=resp,procrustes.p=resp.p)
  return(res)
}

alphaq<-function(mat1,mat2,name.mat1="mat1",name.mat2="mat2",cormethod="pearson",qgrain=1){
  
  
  qz1<-round(seq(from=0, to=5, by=qgrain),digits=1)
  qz2<-round(seq(from=0, to=5, by=qgrain),digits=1)
  
  resm<-matrix(ncol=length(qz1),nrow=length(qz2)) #correlation
  colnames(resm)<-paste(name.mat1,qz2,sep="_")
  rownames(resm)<-paste(name.mat2,qz1,sep="_")
  

  for(i in 1:length(qz1)){
    for(j in 1:length(qz2)){
      a<-hill_taxa(mat1,q=qz1[i])
      b<-hill_taxa(mat2,q=qz2[j])
      
      corr<-cor(a,b,method=cormethod)
     
      
      resm[j,i]<-corr
    }
  }
  res<-list(mantel.cor=resm)
  return(res)
}

pool.spearman<-alphaq(mat1=OTUfungiA,mat2=plantsA)
pool<-alphaq(mat1=OTUfungiA,mat2=plantsA,cormethod="pearson")

plotcorrs<-function(turnqobject){
  obname<-deparse(match.call()$turnqobject)
  coW<-melt(turnqobject$mantel.cor)
  
  ggplot(coW[,-4],aes(Var1,Var2))+ 
    geom_tile(aes(fill = value))+ 
    scale_fill_gradient2(low = "chartreuse4",
                       mid = "white",
                       high = "slateblue4",limits=c(-1,1)) + 
    theme(panel.grid.major.x=element_blank(), 
          panel.grid.minor.x=element_blank(), 
          panel.grid.major.y=element_blank(), 
          panel.grid.minor.y=element_blank(),
          panel.background=element_rect(fill="white"), 
          axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 5,face = "bold"), 
          plot.title = element_text(size=20,face="bold"),
          axis.text.y = element_text(size = 5,face = "bold")) + 
    ggtitle(paste(obname,"(Mantel)",sep="_")) + theme(legend.title=element_text(face="bold", size=14))  + labs(fill="Corr. Coef.")
  }

plotcorrs<-function(turnqobject){
  obname<-deparse(match.call()$turnqobject)
  coW<-melt(turnqobject$mantel.cor)

profileqfunction<-function(datt,qgrain=1){
  qz<-seq(from=0,to=5,by=qgrain)
  res<-matrix(nrow=6,ncol=length(qz))
  for(i in 1:length(qz)){
    temp<-apply(datt,1,d,lev="alpha",q=qz[i])
    res[,i]<-c(qz[i],temp)
  }
  return(res)
}

plotprofiledelta<-function(ambientprofileobject,warmingprofileoject,maintitle="Profile Plot - delta alpha diversity"){
  ppAMFA<-ambientprofileobject
  ppAMFW<-warmingprofileoject
  #the naming of the objects is a legacy for AMF
  topdelta<-max(c(ppAMFA[2,]-ppAMFW[2,],ppAMFA[3,]-ppAMFW[3,],ppAMFA[4,]-ppAMFW[4,],ppAMFA[5,]-ppAMFW[5,],ppAMFA[6,]-ppAMFW[6,]))+5
  bottomdelta<-min(c(ppAMFA[2,]-ppAMFW[2,],ppAMFA[3,]-ppAMFW[3,],ppAMFA[4,]-ppAMFW[4,],ppAMFA[5,]-ppAMFW[5,],ppAMFA[6,]-ppAMFW[6,]))-5
  plot(ppAMFA[1,],ppAMFA[2,],ylim=c(bottomdelta,topdelta),xlim=c(0,6),type="n",ylab="ambient minus warm (effective species)",xlab="Order q",main=maintitle)
  abline(h=0,lty=3)
  lines(ppAMFA[1,],(ppAMFA[2,]-ppAMFW[2,]),col="black");
  lines(ppAMFA[1,],(ppAMFA[3,]-ppAMFW[3,]),col="orange");
  lines(ppAMFA[1,],(ppAMFA[4,]-ppAMFW[4,]),col="red");
  lines(ppAMFA[1,],(ppAMFA[5,]-ppAMFW[5,]),col="blue");
  lines(ppAMFA[1,],(ppAMFA[6,]-ppAMFW[6,]),col="darkgreen");
  
  text(1,topdelta-5,"Plot 1 and 2",cex=0.5,col="black")
  text(2,topdelta-5,"Plot 3 and 4",cex=0.5,col="orange")
  text(3,topdelta-5,"Plot 5 and 6",cex=0.5,col="red")
  text(4,topdelta-5,"Plot 7 and 8",cex=0.5,col="blue")
  text(5,topdelta-5,"Plot 9 and 10",cex=0.5,col="darkgreen")
}




#####


#read in community matrix
plants<-read.csv("~/Desktop/ActiveProjects/WarmingExperiment/CWM/myfileOTUPlants_2023_06_01.csv")[,-1]#delete colum with plot ID
OTUfungi<-t(read.csv("~/Desktop/ActiveProjects/WarmingExperiment/CWM/myfileOTUFungiPlotAvg_2023_06_14.csv",header=TRUE)[,-1])[c(1,3:10,2),];OTUfungi[1:10,1:3] #reorder so 10 is 10th row
OTUAMF<-t(read.csv("~/Desktop/ActiveProjects/WarmingExperiment/CWM/myfileOTUFungiAMFPlotAvg_2023_06_14.csv",header=TRUE)[,-1])[c(1,3:10,2),];OTUAMF[1:10,1:3] #reorder so 10 is 10th row
OTUfungipathogen<-t(read.csv("~/Desktop/ActiveProjects/WarmingExperiment/CWM/myfileOTUFungiPathogensPlotAvg_2023_06_28.csv",header=TRUE)[,-1])[c(1,3:10,2),];OTUfungipathogen[1:10,1:3]#reorder so 10 is 10th row
OTUBacteria<-  t(read.csv("~/Desktop/ActiveProjects/WarmingExperiment/CWM/myfileOTUBacteriaPlotAvg_2023_06_13.csv",header=TRUE)[,-1])[c(1,3:10,2),];OTUBacteria[1:10,1:3]#reorder so 10 is 10th row
OTUSAP<-    t(read.csv("~/Desktop/ActiveProjects/WarmingExperiment/CWM/myfileOTUFungiSaprotrophPlotAvg_2023_06_28.csv",header=TRUE)[,-1])[c(1,3:10,2),];OTUSAP[1:10,1:3]#reorder so 10 is 10th row
OTUDSE<-t(read.csv("~/Desktop/ActiveProjects/WarmingExperiment/CWM/myfileOTUFungiDSEPlotAvg_2023_06_28.csv",header=TRUE)[,-1])[c(1,3:10,2),];OTUDSE[1:10,1:3]#reorder so 10 is 10th row

#partition between warm and ambient
plantsA<-plants[c(1,3,5,7,9),]
plantsW<-plants[c(2,4,6,8,10),]

OTUfungiA<-OTUfungi[c(1,3,5,7,9),]
OTUfungiW<-OTUfungi[c(2,4,6,8,10),]

OTUAMFA<-OTUAMF[c(1,3,5,7,9),]
OTUAMFW<-OTUAMF[c(2,4,6,8,10),]

OTUfungipathogenA<-OTUfungipathogen[c(1,3,5,7,9),]
OTUfungipathogenW<-OTUfungipathogen[c(2,4,6,8,10),]

OTUBacteriaA<-OTUBacteria[c(1,3,5,7,9),]
OTUBacteriaW<-OTUBacteria[c(2,4,6,8,10),]

OTUSAPA<-OTUSAP[c(1,3,5,7,9),]
OTUSAPW<-OTUSAP[c(2,4,6,8,10),]

OTUDSEA<-OTUDSE[c(1,3,5,7,9),]
OTUDSEW<-OTUDSE[c(2,4,6,8,10),]


#questions
#Plants v. AM fungi
 theqgrain<-0.1
plant.AMF.W<-turnq(mat1=plantsW,mat2=OTUAMFW,name.mat1 = "plantsW",name.mat2="AMFW",qgrain=theqgrain)
plant.AMF.A<-turnq(mat1=plantsA,mat2=OTUAMFA,name.mat1 = "plantsA",name.mat2="AMFA",qgrain=theqgrain)
#Plants v. all fungi
 plant.Fungi.W<-turnq(mat1=plantsW,mat2=OTUfungiW,name.mat1 = "plantsW",name.mat2="fungiW",qgrain=theqgrain)
 plant.Fungi.A<-turnq(mat1=plantsA,mat2=OTUfungiA,name.mat1 = "plantsA",name.mat2="fungiA",qgrain=theqgrain)
#Plants v. DSE
 plant.DSE.W<-turnq(mat1=plantsW,mat2=OTUDSEW,name.mat1 = "plantsW",name.mat2="DSEW",qgrain=theqgrain)
 plant.DSE.A<-turnq(mat1=plantsA,mat2=OTUDSEA,name.mat1 = "plantsA",name.mat2="DSEA",qgrain=theqgrain)
 #Plants v. Pathogens
 plant.path.W<-turnq(mat1=plantsW,mat2=OTUfungipathogenW,name.mat1 = "plantsW",name.mat2="PathW",qgrain=theqgrain)
 plant.path.A<-turnq(mat1=plantsA,mat2=OTUfungipathogenA,name.mat1 = "plantsA",name.mat2="PathA",qgrain=theqgrain)
 #Plants v. Saprotrophs
 plant.SAP.W<-turnq(mat1=plantsW,mat2=OTUSAPW,name.mat1 = "plantsW",name.mat2="SAPW",qgrain=theqgrain)
 plant.SAP.A<-turnq(mat1=plantsA,mat2=OTUSAPA,name.mat1 = "plantsA",name.mat2="SAPA",qgrain=theqgrain)
 #AM fungi v. DSE
 AMF.DSE.W<-turnq(mat1=OTUAMFW,mat2=OTUDSEW,name.mat1 = "AMFW",name.mat2="DSEW",qgrain=theqgrain)
 AMF.DSE.A<-turnq(mat1=OTUAMFA,mat2=OTUDSEA,name.mat1 = "AMFA",name.mat2="DSEA",qgrain=theqgrain)
 
#AM fungi v. Saprotroph
 AMF.SAP.W<-turnq(mat1=OTUAMFW,mat2=OTUSAPW,name.mat1 = "AMFW",name.mat2="SAPW",qgrain=theqgrain)
 AMF.SAP.A<-turnq(mat1=OTUAMFA,mat2=OTUSAPA,name.mat1 = "AMFA",name.mat2="SAPA",qgrain=theqgrain)



mantel.correlations<-list(plant.AMF.W=plant.AMF.W,plant.AMF.A=plant.AMF.A,plant.Fungi.W=plant.Fungi.W,plant.Fungi.A=plant.Fungi.A,plant.DSE.W=plant.DSE.W,plant.DSE.A=plant.DSE.A,plant.path.W=plant.path.W,plant.path.A=plant.path.A,plant.SAP.W=plant.SAP.W,plant.SAP.A=plant.SAP.A,AMF.DSE.W=AMF.DSE.W,AMF.DSE.A=AMF.DSE.A,AMF.SAP.W=AMF.SAP.W,AMF.SAP.A=AMF.SAP.A)
AMF.SAP.A$mantel.cor
load("mantel.correlations.Rdata") #above code (~30 minutes in terminal)
plant.AMF.W<-mantel.correlations$plant.AMF.W
plant.AMF.A<-mantel.correlations$plant.AMF.A
plant.Fungi.W<-mantel.correlations$plant.Fungi.W
plant.Fungi.A<-mantel.correlations$plant.Fungi.A
plant.DSE.W<-mantel.correlations$plant.DSE.W
plant.DSE.A<-mantel.correlations$plant.DSE.A
plant.path.W<-mantel.correlations$plant.path.W
plant.path.A<-mantel.correlations$plant.path.A
plant.SAP.W<-mantel.correlations$plant.SAP.W
plant.SAP.A<-mantel.correlations$plant.SAP.A
AMF.DSE.W<-mantel.correlations$AMF.DSE.W
AMF.DSE.A<-mantel.correlations$AMF.DSE.A
AMF.SAP.W<-mantel.correlations$AMF.SAP.W
AMF.SAP.A<-mantel.correlations$AMF.SAP.A

#plot them
plotcorrs(plant.AMF.W)
plotcorrs(plant.AMF.A)

plotcorrs(plant.Fungi.W)
plotcorrs(plant.Fungi.A)

plotcorrs(plant.DSE.W)
plotcorrs(plant.DSE.A)

plotcorrs(plant.path.W)
plotcorrs(plant.path.A)

plotcorrs(plant.SAP.W)
plotcorrs(plant.SAP.A)

plotcorrs(AMF.DSE.W)
plotcorrs(AMF.DSE.A)

plotcorrs(AMF.SAP.W)
plotcorrs(AMF.SAP.A)

#Alpha Diversity
#plot delta alpha  Ambient minus Warming. So, positive numbers indicate more alpha diversity in ambient than warming.

pp.plant.W<-profileqfunction(plantsW,qgrain=0.1)
pp.plant.A<-profileqfunction(plantsA,qgrain=0.1)
plotprofiledelta(pp.plant.A,pp.plant.W,maintitle = "Plant - delta alpha")

pp.AMF.W<-profileqfunction(OTUAMFW,qgrain=0.1)
pp.AMF.A<-profileqfunction(OTUAMFA,qgrain=0.1)
plotprofiledelta(pp.AMF.A,pp.AMF.W,maintitle = "AMF - delta alpha")

pp.fungi.W<-profileqfunction(OTUfungiW,qgrain=0.1)
pp.fungi.A<-profileqfunction(OTUfungiA,qgrain=0.1)
plotprofiledelta(pp.fungi.A,pp.fungi.W,maintitle = "Fungi all - delta alpha")

pp.DSE.W<-profileqfunction(OTUDSEW,qgrain=0.1)
pp.DSE.A<-profileqfunction(OTUDSEA,qgrain=0.1)
plotprofiledelta(pp.DSE.A,pp.DSE.W,maintitle = "DSE - delta alpha")

pp.pathogen.W<-profileqfunction(OTUfungipathogenW,qgrain=0.1)
pp.pathogen.A<-profileqfunction(OTUfungipathogenA,qgrain=0.1)
plotprofiledelta(pp.pathogen.A,pp.pathogen.W,maintitle = "Pathogen - delta alpha")

pp.SAP.W<-profileqfunction(OTUSAPW,qgrain=0.1)
pp.SAP.A<-profileqfunction(OTUSAPA,qgrain=0.1)
plotprofiledelta(pp.SAP.A,pp.SAP.W,maintitle = "SAP - delta alpha")


pp.bacteria.W<-profileqfunction(OTUBacteriaW,qgrain=0.1)
pp.bacteria.A<-profileqfunction(OTUBacteriaA,qgrain=0.1)
plotprofiledelta(pp.bacteria.A,pp.bacteria.W,maintitle = "Bacteria - delta alpha")

theqgrain<-0.1
plant.AMF.W<-alphaq(mat1=plantsW,mat2=OTUAMFW,name.mat1 = "plantsW",name.mat2="AMFW",qgrain=theqgrain)
plant.AMF.A<-alphaq(mat1=plantsA,mat2=OTUAMFA,name.mat1 = "plantsA",name.mat2="AMFA",qgrain=theqgrain)
#Plants v. all fungi
 plant.Fungi.W<-alphaq(mat1=plantsW,mat2=OTUfungiW,name.mat1 = "plantsW",name.mat2="fungiW",qgrain=theqgrain)
 plant.Fungi.A<-alphaq(mat1=plantsA,mat2=OTUfungiA,name.mat1 = "plantsA",name.mat2="fungiA",qgrain=theqgrain)
#Plants v. DSE
 plant.DSE.W<-alphaq(mat1=plantsW,mat2=OTUDSEW,name.mat1 = "plantsW",name.mat2="DSEW",qgrain=theqgrain)
 plant.DSE.A<-alphaq(mat1=plantsA,mat2=OTUDSEA,name.mat1 = "plantsA",name.mat2="DSEA",qgrain=theqgrain)
 #Plants v. Pathogens
 plant.path.W<-alphaq(mat1=plantsW,mat2=OTUfungipathogenW,name.mat1 = "plantsW",name.mat2="PathW",qgrain=theqgrain)
 plant.path.A<-alphaq(mat1=plantsA,mat2=OTUfungipathogenA,name.mat1 = "plantsA",name.mat2="PathA",qgrain=theqgrain)
 #Plants v. Saprotrophs
 plant.SAP.W<-alphaq(mat1=plantsW,mat2=OTUSAPW,name.mat1 = "plantsW",name.mat2="SAPW",qgrain=theqgrain)
 plant.SAP.A<-alphaq(mat1=plantsA,mat2=OTUSAPA,name.mat1 = "plantsA",name.mat2="SAPA",qgrain=theqgrain)
 #AM fungi v. DSE
 AMF.DSE.W<-alphaq(mat1=OTUAMFW,mat2=OTUDSEW,name.mat1 = "AMFW",name.mat2="DSEW",qgrain=theqgrain)
 AMF.DSE.A<-alphaq(mat1=OTUAMFA,mat2=OTUDSEA,name.mat1 = "AMFA",name.mat2="DSEA",qgrain=theqgrain)
 
#AM fungi v. Saprotroph
 AMF.SAP.W<-alphaq(mat1=OTUAMFW,mat2=OTUSAPW,name.mat1 = "AMFW",name.mat2="SAPW",qgrain=theqgrain)
 AMF.SAP.A<-alphaq(mat1=OTUAMFA,mat2=OTUSAPA,name.mat1 = "AMFA",name.mat2="SAPA",qgrain=theqgrain)

##### ALPHA DIVERSITY STUFF NOTES #####
library(hillR);library(vegan);library(vegetarian)

#read in community matrix
plants<-read.csv("myfileOTUPlants_2023_06_01.csv")[,-1]#delete colum with plot ID
OTUfungi<-t(read.csv("myfileOTUFungiPlotAvg_2023_06_14.csv",header=TRUE)[,-1])[c(1,3:10,2),] #reorder so 10 is 10th row
OTUAMF<-t(read.csv("myfileOTUFungiAMFPlotAvg_2023_06_14.csv",header=TRUE)[,-1])[c(1,3:10,2),] #reorder so 10 is 10th row
OTUfungipathogen<-t(read.csv("myfileOTUFungiPathogensPlotAvg_2023_06_14.csv",header=TRUE)[,-1])[c(1,3:10,2),]#reorder so 10 is 10th row
OTUBacteria<-  t(read.csv("myfileOTUBacteriaPlotAvg_2023_06_13.csv",header=TRUE)[,-1])[c(1,3:10,2),]#reorder so 10 is 10th row
OTUSAP<-    t(read.csv("myfileOTUFungiSaprotrophsPlotAvg_2023_06_14.csv",header=TRUE)[,-1])[c(1,3:10,2),]

#partition between warm and ambient
plantsA<-plants[c(1,3,5,7,9),]
plantsW<-plants[c(2,4,6,8,10),]
OTUfungiA<-OTUfungi[c(1,3,5,7,9),]
OTUfungiW<-OTUfungi[c(2,4,6,8,10),]
OTUAMFA<-OTUAMF[c(1,3,5,7,9),]
OTUAMFW<-OTUAMF[c(2,4,6,8,10),]
OTUfungipathogenA<-OTUfungipathogen[c(1,3,5,7,9),]
OTUfungipathogenW<-OTUfungipathogen[c(2,4,6,8,10),]
OTUBacteriaA<-OTUBacteria[c(1,3,5,7,9),]
OTUBacteriaW<-OTUBacteria[c(2,4,6,8,10),]
OTUSAPA<-OTUSAP[c(1,3,5,7,9),]
OTUSAPW<-OTUSAP[c(2,4,6,8,10),]


ppAMFA<-profileqfunction(OTUAMFA,qgrain=0.1)
ppAMFW<-profileqfunction(OTUAMFW,qgrain=0.1)


top<-max(rbind(ppAMFA,ppAMFW))+10
#plot individual plots
plot(ppAMFA[1,],ppAMFA[2,],ylim=c(0,top),type="n")
lines(ppAMFA[1,],ppAMFA[2,])
lines(ppAMFA[1,],ppAMFA[3,])
lines(ppAMFA[1,],ppAMFA[4,])
lines(ppAMFA[1,],ppAMFA[5,])
lines(ppAMFA[1,],ppAMFA[6,])

lines(ppAMFA[1,],ppAMFW[2,],col="red")
lines(ppAMFA[1,],ppAMFW[3,],col="red")
lines(ppAMFA[1,],ppAMFW[4,],col="red")
lines(ppAMFA[1,],ppAMFW[5,],col="red")
lines(ppAMFA[1,],ppAMFW[6,],col="red")


















