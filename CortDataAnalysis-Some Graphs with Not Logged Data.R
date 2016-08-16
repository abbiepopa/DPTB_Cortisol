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
AllDataCortAnalysis$CortRatio<-AllDataCortAnalysis$PostCORT/AllDataCortAnalysis$PreCORT

###Plots to understand cort values better
quartz()
ggplot(AllDataCortAnalysis, aes(factor(Dx),PreCORT))+geom_violin(trim = F, aes(fill=factor(Dx)))+scale_fill_manual(values=c("maroon2","blue"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Diagnosis") + ylab("Pre Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(Dx),PostCORT))+geom_violin(trim = F, aes(fill=factor(Dx)))+scale_fill_manual(values=c("maroon2","blue"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Diagnosis") + ylab("Post Cortisol")

###scatters comparing delta, ratio, and log-log
quartz()
plot(AllDataCortAnalysis$CortLogDelta, AllDataCortAnalysis$CortDelta, col=ifelse(AllDataCortAnalysis$Dx=="22q","maroon2","blue"), ylab="Post-Pre", xlab="LogPost-LogPre")

plot(AllDataCortAnalysis$CortRatio, AllDataCortAnalysis$CortDelta, col=ifelse(AllDataCortAnalysis$Dx=="22q","maroon2","blue"), ylab="Post-Pre", xlab="Post/Pre")

###significant scatters but generate the scatters on unlogged###
quartz()
plot(AllDataCortAnalysis$PreCORT, AllDataCortAnalysis$PercHappy, main="r=0.32 p=0.04", xlab="PreCort", ylab="PercHappy",col=ifelse(AllDataCortAnalysis$Dx=="22q","maroon2","blue"))

quartz()
plot(AllDataCortAnalysis$PreCORT, AllDataCortAnalysis$PercHappyNonFAce, main="r=-0.32 p=0.05", xlab="PreCort", ylab="PercHappyNonFace",col=ifelse(AllDataCortAnalysis$Dx=="22q","maroon2","blue"))

###honestly this severely seems to be not working. It might be useful to look at spaghetti without logs though

CortDxWide<-AllDataCortAnalysis[,c(1,2,8,9,77)]
library(reshape)
CortDxLong<-melt(CortDxWide, id=c("CABIL_ID","Dx", "CortDelta"))
colnames(CortDxLong)[4]<-"Time"

CortDxLongOrdered<-CortDxLong[order(CortDxLong$CABIL_ID),]
CortDxLongOrderedNoNA<-na.omit(CortDxLongOrdered)

quartz()
interaction.plot(CortDxLongOrderedNoNA$Time,trace.factor=CortDxLongOrderedNoNA$CABIL_ID, response=CortDxLongOrderedNoNA$value, legend=F, xlab="Time", ylab="Log Cortisol", ylim=c(0,2.5), xaxt=F)
+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))

LongCort22q<-subset(CortDxLongOrderedNoNA, Dx=="22q")
LongCortTD<-subset(CortDxLongOrderedNoNA, Dx=="TD")

quartz()
interaction.plot(LongCort22q$Time,trace.factor=LongCort22q$CABIL_ID, response=LongCort22q$value, legend=F, col="maroon2", xlab="Time", ylab="Log Cortisol", ylim=c(0,2.5), xaxt=F)
+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))
quartz()
interaction.plot(LongCortTD$Time,trace.factor=LongCortTD$CABIL_ID, response=LongCortTD$value, legend=F, col="blue", xlab="Time", ylab="Log Cortisol", ylim=c(0,2.5), xaxt=F)
+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))