setwd("/Users/abbiepopa/Documents/Lab/DPTB/Cortisol Analysis/cortdata from Elliot/Working")
Participants<-read.csv("PartList_DPTB.csv", na.strings=".")

RTTD<-read.csv("DPTB_RT_9-17-14_TDKoralyOut.csv",na.strings="NaN")
RTTD<-RTTD[,c(1,9,10)]

RT22q<-read.csv("DPTB_RT_9-17-14_22qKoralyOut.csv", na.strings="NaN")
RT22q<-RT22q[,c(1,9,10)]

RT<-rbind(RTTD,RT22q)


Cort<-read.csv("CORT-Data-EAB_July2014.csv")
Cort<-Cort[,c(1, 3:6)]
colnames(Cort)[1]<-"CABIL_ID"

SpenceABAS<-read.csv("Clustering_Database_8-8-14.csv", na.strings="-9999")

EyeGazeOverall<-read.csv("AllData.csv", na.strings="NA")
EyeGazeOverall<-EyeGazeOverall[,2:11]
colnames(EyeGazeOverall)[1]<-"CABIL_ID"

EyeGazeTimeCourse<-read.csv("AllDataTC25no676.csv", na.strings="NA")
EyeGazeTimeCourse<-EyeGazeTimeCourse[,2:24]
colnames(EyeGazeTimeCourse)[1]<-"CABIL_ID"

PupilOverall22q<-read.csv("PupilPilot_6-5-14_22q.csv", na.strings=".")
PupilOverall22q<-PupilOverall22q[,c(1:7, 9:13)]
colnames(PupilOverall22q)[1]<-"CABIL_ID"

PupilOverallTD<-read.csv("PupilPilot_6-5-14_TD.csv", na.strings=".")
PupilOverallTD<-PupilOverallTD[,c(1:7,9:13)]
colnames(PupilOverallTD)[1]<-"CABIL_ID"

PupilOverall<-rbind(PupilOverall22q, PupilOverallTD)

PupilChangeTD<-read.csv("PupilChange_6-5-14_TD.csv",na.strings=".")
PupilChangeTD<-PupilChangeTD[,1:7]
colnames(PupilChangeTD)[1]<-"CABIL_ID"

PupilChange22q<-read.csv("PupilChange_6-5-14_22q.csv",na.strings=".")
PupilChange22q<-PupilChange22q[,1:7]
colnames(PupilChange22q)[1]<-"CABIL_ID"

PupilChange<-rbind(PupilChangeTD,PupilChange22q)

AllDataCortAnalysis<-merge(Participants, Cort, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, SpenceABAS, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, RT, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, EyeGazeOverall, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, EyeGazeTimeCourse, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, PupilOverall, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, PupilChange, all.x=T)

AllDataCortAnalysis$CortDelta<-AllDataCortAnalysis$PostCORT-AllDataCortAnalysis$PreCORT
AllDataCortAnalysis$LogCortDelta<-log10(AllDataCortAnalysis$CortDelta+1)
AllDataCortAnalysis$CortLogDelta<-AllDataCortAnalysis$LogPostCort-AllDataCortAnalysis$LogPreCort

###clusters###
ClusAll<-read.csv("AllDataClus.csv")
Clus<-ClusAll[,c(2,12)]
colnames(Clus)[1]<-"CABIL_ID"

AllDataCortAnalysis<-merge(AllDataCortAnalysis, Clus, all.x=T)

###Correlations, everyone together

for(i in 5:76){
	pretest<-cor.test(AllDataCortAnalysis$LogPreCort, AllDataCortAnalysis[,i])
	if(pretest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Cluster==2,"green","purple"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(AllDataCortAnalysis$LogPostCort, AllDataCortAnalysis[,i])
	if(posttest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPostCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Cluster==2,"green","purple"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(AllDataCortAnalysis$CortDelta, AllDataCortAnalysis[,i])
	if(deltatest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$Cluster==2,"green","purple"))
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

###Cluster 1

ClusterOne<-subset(AllDataCortAnalysis, Cluster==1)

for(i in 5:76){
	pretest<-cor.test(ClusterOne$LogPreCort, ClusterOne[,i])
	if(pretest$p.value<0.05){
			nowData<-ClusterOne[which(!is.na(ClusterOne$LogPreCort)&!is.na(ClusterOne[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Cluster==2,"green","purple"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(ClusterOne$LogPostCort, ClusterOne[,i])
	if(posttest$p.value<0.05){
			nowData<-ClusterOne[which(!is.na(ClusterOne$LogPostCort)&!is.na(ClusterOne[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Cluster==2,"green","purple"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(ClusterOne$CortDelta, ClusterOne[,i])
	if(deltatest$p.value<0.05){
			nowData<-ClusterOne[which(!is.na(ClusterOne$CortDelta)&!is.na(ClusterOne[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$Cluster==2,"green","purple"))
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

###Cluster 2

ClusterTwo<-subset(AllDataCortAnalysis, Cluster==2)

for(i in 5:76){
	pretest<-cor.test(ClusterTwo$LogPreCort, ClusterTwo[,i])
	if(pretest$p.value<0.05){
			nowData<-ClusterTwo[which(!is.na(ClusterTwo$LogPreCort)&!is.na(ClusterTwo[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Cluster==2,"green","purple"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(ClusterTwo$LogPostCort, ClusterTwo[,i])
	if(posttest$p.value<0.05){
			nowData<-ClusterTwo[which(!is.na(ClusterTwo$LogPostCort)&!is.na(ClusterTwo[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Cluster==2,"green","purple"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(ClusterTwo$CortDelta, ClusterTwo[,i])
	if(deltatest$p.value<0.05){
			nowData<-ClusterTwo[which(!is.na(ClusterTwo$CortDelta)&!is.na(ClusterTwo[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$Cluster==2,"green","purple"))
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

###group t-tests
#plot(AllDataCortAnalysis$Clus, AllDataCortAnalysis$LogPreCort)
t.test(AllDataCortAnalysis$LogPreCort[which(AllDataCortAnalysis$Cluster==1)], AllDataCortAnalysis$LogPreCort[which(AllDataCortAnalysis$Cluster==2)])
t.test(AllDataCortAnalysis$LogPostCort[which(AllDataCortAnalysis$Cluster==1)], AllDataCortAnalysis$LogPostCort[which(AllDataCortAnalysis$Cluster==2)])
t.test(AllDataCortAnalysis$CortDelta[which(AllDataCortAnalysis$Cluster==1)], AllDataCortAnalysis$CortDelta[which(AllDataCortAnalysis$Cluster==2)])
t.test(AllDataCortAnalysis$LogCortDelta[which(AllDataCortAnalysis$Cluster==1)], AllDataCortAnalysis$LogCortDelta[which(AllDataCortAnalysis$Cluster==2)])

library(ggplot2)
###Plots to understand cort values better
ClusteredIDs<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$Cluster)),]
ClusteredIDs$ClusName<-ifelse(ClusteredIDs$Cluster==1, "Struggler", "Coper")
quartz()
ggplot(ClusteredIDs, aes(factor(ClusName),LogPreCort))+geom_violin(trim = F, aes(fill=factor(ClusName)))+scale_fill_manual(values=c("green","purple"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Cluster") + ylab("Log Pre Cortisol")
quartz()
ggplot(ClusteredIDs, aes(factor(ClusName),LogPostCort))+geom_violin(trim = F, aes(fill=factor(ClusName)))+scale_fill_manual(values=c("green","purple"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Cluster") + ylab("Log Post Cortisol")
quartz()
ggplot(ClusteredIDs, aes(factor(ClusName),CortDelta))+geom_violin(trim = F, aes(fill=factor(ClusName)))+scale_fill_manual(values=c("green","purple"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Cluster") + ylab("Change in Cortisol")

quartz()
ggplot(ClusteredIDs, aes(factor(ClusName),PreCORT))+geom_violin(trim = F, aes(fill=factor(ClusName)))+scale_fill_manual(values=c("green","purple"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Cluster") + ylab("Pre Cortisol")
quartz()
ggplot(ClusteredIDs, aes(factor(ClusName),PostCORT))+geom_violin(trim = F, aes(fill=factor(ClusName)))+scale_fill_manual(values=c("green","purple"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Cluster") + ylab("Post Cortisol")

CortClusWideLog<-ClusteredIDs[,c(1,10,11,81)]
CortClusWide<-ClusteredIDs[,c(1,8,9,81)]
library(reshape)
CortClusLongLog<-melt(CortClusWideLog, id=c("CABIL_ID","ClusName"))
CortClusLong<-melt(CortClusWide, id=c("CABIL_ID","ClusName"))

colnames(CortClusLongLog)[3]<-"Time"
colnames(CortClusLong)[3]<-"Time"

LogLongCortCoper<-subset(CortClusLongLog, ClusName=="Coper")
LogLongCortStruggler<-subset(CortClusLongLog, ClusName=="Struggler")

LongCortCoper<-subset(CortClusLong, ClusName=="Coper")
LongCortStruggler<-subset(CortClusLong, ClusName=="Struggler")

naoLogLongCortCoper<-na.omit(LogLongCortCoper)
naoLogLongCortStruggler<-na.omit(LogLongCortStruggler)

snaoLogLongCortCoper<-naoLogLongCortCoper[order(naoLogLongCortCoper$CABIL_ID),]

quartz()
interaction.plot(snaoLogLongCortCoper$Time,trace.factor=snaoLogLongCortCoper$CABIL_ID, response=snaoLogLongCortCoper$value, legend=F, col="green", xlab="Time", ylab="Log Cortisol", ylim=c(-1.5,0.5), xaxt=F)
+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))
quartz()
interaction.plot(naoLogLongCortStruggler$Time,trace.factor=naoLogLongCortStruggler$CABIL_ID, response=naoLogLongCortStruggler$value, legend=F, col="purple", xlab="Time", ylab="Log Cortisol", ylim=c(-1.5,0.5), xaxt=F)
+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))
