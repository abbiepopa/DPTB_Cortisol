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

###note, use log10 

###Correlations, everyone together

for(i in 5:76){
	pretest<-cor.test(AllDataCortAnalysis$LogPreCort, AllDataCortAnalysis[,i])
	if(pretest$p.value<0.05){
	quartz()
	plot(AllDataCortAnalysis$LogPreCort, AllDataCortAnalysis[,i], main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2)), xlab="LogPreCort", ylab=colnames(AllDataCortAnalysis)[i],col=ifelse(AllDataCortAnalysis$Dx=="22q","maroon2","blue"))
	#abline(lm(AllDataCortAnalysis$LogPreCort~AllDataCortAnalysis[,i]))
	}
	posttest<-cor.test(AllDataCortAnalysis$LogPostCort, AllDataCortAnalysis[,i])
	if(posttest$p.value<0.05){
	quartz()
	plot(AllDataCortAnalysis$LogPostCort, AllDataCortAnalysis[,i], main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2)), xlab="LogPostCort", ylab=colnames(AllDataCortAnalysis)[i],col=ifelse(AllDataCortAnalysis$Dx=="22q","maroon2","blue"))
	#abline(lm(AllDataCortAnalysis$LogPostCort~AllDataCortAnalysis[,i]))
	}
	deltatest<-cor.test(AllDataCortAnalysis$CortLogDelta, AllDataCortAnalysis[,i])
	if(deltatest$p.value<0.05){
	quartz()
	plot(AllDataCortAnalysis$CortLogDelta, AllDataCortAnalysis[,i], main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2)), xlab="Change in Cort (Post-Pre)", ylab=colnames(AllDataCortAnalysis)[i],col=ifelse(AllDataCortAnalysis$Dx=="22q","maroon2","blue"))
	#abline(lm(AllDataCortAnalysis$CortLogDelta~AllDataCortAnalysis[,i]))
	}
	# logdeltatest<-cor.test(AllDataCortAnalysis$LogCortLogDelta, AllDataCortAnalysis[,i])
	# if(logdeltatest$p.value<0.05){
	# quartz()
	# plot(AllDataCortAnalysis$LogCortLogDelta, AllDataCortAnalysis[,i], main=paste("r=",round(logdeltatest$estimate, digits=2),"p=",round(logdeltatest$p.value, digits=2)), xlab="Log: Change in Cort (Post-Pre)", ylab=colnames(AllDataCortAnalysis)[i],col=ifelse(AllDataCortAnalysis$Dx=="22q","maroon2","blue"))
	# #abline(lm(AllDataCortAnalysis$LogCortLogDelta~AllDataCortAnalysis[,i]))
	# }
}

###change t-test, everone together
t.test(AllDataCortAnalysis$CortLogDelta)
t.test(AllDataCortAnalysis$LogCortLogDelta)

###group t-tests
plot(AllDataCortAnalysis$Dx, AllDataCortAnalysis$LogPreCort)
t.test(AllDataCortAnalysis$LogPreCort[which(AllDataCortAnalysis$Dx=="22q")], AllDataCortAnalysis$LogPreCort[which(AllDataCortAnalysis$Dx=="TD")])
t.test(AllDataCortAnalysis$LogPostCort[which(AllDataCortAnalysis$Dx=="22q")], AllDataCortAnalysis$LogPostCort[which(AllDataCortAnalysis$Dx=="TD")])
t.test(AllDataCortAnalysis$CortLogDelta[which(AllDataCortAnalysis$Dx=="22q")], AllDataCortAnalysis$CortLogDelta[which(AllDataCortAnalysis$Dx=="TD")])
t.test(AllDataCortAnalysis$LogCortLogDelta[which(AllDataCortAnalysis$Dx=="22q")], AllDataCortAnalysis$LogCortLogDelta[which(AllDataCortAnalysis$Dx=="TD")])

###would change in looking times be a better measure than early and late looking times?

###Correlations, 22q
Cort22q<-subset(AllDataCortAnalysis, Dx=="22q")

for(i in 5:76){
	pretest<-cor.test(Cort22q$LogPreCort, Cort22q[,i])
	if(pretest$p.value<0.05){
	quartz()
	plot(Cort22q$LogPreCort, Cort22q[,i], main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2)), xlab="LogPreCort", ylab=colnames(Cort22q)[i],col=ifelse(Cort22q$Dx=="22q","maroon2","blue"))
	#abline(lm(Cort22q$LogPreCort~Cort22q[,i]))
	}
	posttest<-cor.test(Cort22q$LogPostCort, Cort22q[,i])
	if(posttest$p.value<0.05){
	quartz()
	plot(Cort22q$LogPostCort, Cort22q[,i], main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2)), xlab="LogPostCort", ylab=colnames(Cort22q)[i],col=ifelse(Cort22q$Dx=="22q","maroon2","blue"))
	#abline(lm(Cort22q$LogPostCort~Cort22q[,i]))
	}
	deltatest<-cor.test(Cort22q$CortLogDelta, Cort22q[,i])
	if(deltatest$p.value<0.05){
	quartz()
	plot(Cort22q$CortLogDelta, Cort22q[,i], main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2)), xlab="Change in Cort (Post-Pre)", ylab=colnames(Cort22q)[i],col=ifelse(Cort22q$Dx=="22q","maroon2","blue"))
	#abline(lm(Cort22q$CortLogDelta~Cort22q[,i]))
	}
	logdeltatest<-cor.test(Cort22q$LogCortLogDelta, Cort22q[,i])
	if(logdeltatest$p.value<0.05){
	quartz()
	plot(Cort22q$LogCortLogDelta, Cort22q[,i], main=paste("r=",round(logdeltatest$estimate, digits=2),"p=",round(logdeltatest$p.value, digits=2)), xlab="Log: Change in Cort (Post-Pre)", ylab=colnames(Cort22q)[i],col=ifelse(Cort22q$Dx=="22q","maroon2","blue"))
	#abline(lm(Cort22q$LogCortLogDelta~Cort22q[,i]))
	}
}

###change t-test, 22q
t.test(Cort22q$CortLogDelta)
t.test(Cort22q$LogCortLogDelta)


###Correlations, TD
CortTD<-subset(AllDataCortAnalysis, Dx=="TD")

for(i in 5:76){
	pretest<-cor.test(CortTD$LogPreCort, CortTD[,i])
	if(pretest$p.value<0.05){
	quartz()
	plot(CortTD$LogPreCort, CortTD[,i], main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2)), xlab="LogPreCort", ylab=colnames(CortTD)[i],col=ifelse(CortTD$Dx=="22q","maroon2","blue"))
	#abline(lm(CortTD$LogPreCort~CortTD[,i]))
	}
	posttest<-cor.test(CortTD$LogPostCort, CortTD[,i])
	if(posttest$p.value<0.05){
	quartz()
	plot(CortTD$LogPostCort, CortTD[,i], main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2)), xlab="LogPostCort", ylab=colnames(CortTD)[i],col=ifelse(CortTD$Dx=="22q","maroon2","blue"))
	#abline(lm(CortTD$LogPostCort~CortTD[,i]))
	}
	deltatest<-cor.test(CortTD$CortLogDelta, CortTD[,i])
	if(deltatest$p.value<0.05){
	quartz()
	plot(CortTD$CortLogDelta, CortTD[,i], main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2)), xlab="Change in Cort (Post-Pre)", ylab=colnames(CortTD)[i],col=ifelse(CortTD$Dx=="22q","maroon2","blue"))
	#abline(lm(CortTD$CortLogDelta~CortTD[,i]))
	}
	logdeltatest<-cor.test(CortTD$LogCortLogDelta, CortTD[,i])
	if(logdeltatest$p.value<0.05){
	quartz()
	plot(CortTD$LogCortLogDelta, CortTD[,i], main=paste("r=",round(logdeltatest$estimate, digits=2),"p=",round(logdeltatest$p.value, digits=2)), xlab="Log: Change in Cort (Post-Pre)", ylab=colnames(CortTD)[i],col=ifelse(CortTD$Dx=="22q","maroon2","blue"))
	#abline(lm(CortTD$LogCortLogDelta~CortTD[,i]))
	}
}

###change t-test, TD
t.test(CortTD$CortLogDelta)
t.test(CortTD$LogCortLogDelta)


###Plots to understand cort values better
quartz()
ggplot(AllDataCortAnalysis, aes(factor(Dx),LogPreCort))+geom_violin(trim = F, aes(fill=factor(Dx)))+scale_fill_manual(values=c("maroon2","blue"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Diagnosis") + ylab("Log Pre Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(Dx),LogPostCort))+geom_violin(trim = F, aes(fill=factor(Dx)))+scale_fill_manual(values=c("maroon2","blue"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Diagnosis") + ylab("Log Post Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(Dx),CortLogDelta))+geom_violin(trim = F, aes(fill=factor(Dx)))+scale_fill_manual(values=c("maroon2","blue"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Diagnosis") + ylab("Change in Cortisol")

CortDxWide<-AllDataCortAnalysis[,c(1,2,10,11,77)]
library(reshape)
CortDxLong<-melt(CortDxWide, id=c("CABIL_ID","Dx", "CortLogDelta"))
colnames(CortDxLong)[4]<-"Time"

quartz()
interaction.plot(CortDxLong$Time,trace.factor=CortDxLong$CABIL_ID, response=CortDxLong$value, legend=F, col=ifelse(CortDxLong$Dx=="22q","maroon2","blue"))

CortDxLongOrdered<-CortDxLong[order(CortDxLong$CABIL_ID),]
# CortDxLongOrderedNoNA<-na.omit(CortDxLongOrdered)

# DxColoring<-function(x){
	# if(x=="22q"){
		# mycolor<-"maroon2"} else{if(x=="TD"){
			# mycolor<-"blue"} else {
				# mycolor<-"gray"}}
# return(mycolor)}
# MyDxColorVector<-matrix(NA, nrow=170, ncol=1)
# for(i in 1:170){
	# MyDxColorVector[i]<-DxColoring(CortDxLongOrderedNoNA$Dx[i])
# }

# #MyDxColorVector<-MyDxColorVector[seq(from=1, to=169, by=2),]

# RxColoring<-function(x){
	# if(is.na(x)){
		# mycolor<-"blue"
	# }else{
	# if(abs(x)>0.25){
		# mycolor<-"red"} else {
				# mycolor<-"gray"}}
# return(mycolor)}
# MyRxColorVector<-matrix(NA, nrow=170, ncol=1)
# for(i in 1:170){
	# MyRxColorVector[i]<-RxColoring(CortDxLongOrderedNoNA$CortLogDelta[i])
	# i
# }

# #MyRxColorVector<-MyRxColorVector[seq(from=1, to=169, by=2),]

# quartz()
# interaction.plot(CortDxLongOrderedNoNA$Time,trace.factor=CortDxLongOrderedNoNA$CABIL_ID, response=CortDxLongOrderedNoNA$value, legend=T, col=MyRxColorVector)

# quartz()
# interaction.plot(CortDxLongOrderedNoNA$Time,trace.factor=CortDxLongOrderedNoNA$CABIL_ID, response=CortDxLongOrderedNoNA$value, legend=T, col=MyDxColorVector)

quartz()
interaction.plot(CortDxLongOrderedNoNA$Time,trace.factor=CortDxLongOrderedNoNA$CABIL_ID, response=CortDxLongOrderedNoNA$value, legend=F, xlab="Time", ylab="Log Cortisol", ylim=c(-1.5,0.5), xaxt=F)
+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))

LongCort22q<-subset(CortDxLongOrderedNoNA, Dx=="22q")
LongCortTD<-subset(CortDxLongOrderedNoNA, Dx=="TD")

CortLogDeltaMean<-mean(CortDxLongOrderedNoNA$CortLogDelta)
CortLogDeltaSD<-sd(CortDxLongOrderedNoNA$CortLogDelta)
CortLogDelta1Up<-CortLogDeltaMean+CortLogDeltaSD
CortLogDelta1Down<-CortLogDeltaMean-CortLogDeltaSD
CortLogDelta2Up<-CortLogDeltaMean+(2*CortLogDeltaSD)
CortLogDelta2Down<-CortLogDeltaMean-(2*CortLogDeltaSD)

LongCortBigBigRx<-subset(CortDxLongOrderedNoNA, CortLogDelta>CortLogDelta2Up | CortLogDelta<CortLogDelta2Down)
LongCortBigRx<-subset(CortDxLongOrderedNoNA, CortLogDelta>CortLogDelta1Up & CortLogDelta<CortLogDelta2Up | CortLogDelta<CortLogDelta1Down & CortLogDelta>CortLogDelta2Down)
LongCortLittleRx<-subset(CortDxLongOrderedNoNA, CortLogDelta<CortLogDelta1Up & CortLogDelta>CortLogDelta1Down)

quartz()
interaction.plot(LongCort22q$Time,trace.factor=LongCort22q$CABIL_ID, response=LongCort22q$value, legend=F, col="maroon2", xlab="Time", ylab="Log Cortisol", ylim=c(-1.5,0.5), xaxt=F)
+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))
quartz()
interaction.plot(LongCortTD$Time,trace.factor=LongCortTD$CABIL_ID, response=LongCortTD$value, legend=F, col="blue", xlab="Time", ylab="Log Cortisol", ylim=c(-1.5,0.5), xaxt=F)
+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))
quartz()
interaction.plot(LongCortBigBigRx$Time,trace.factor=LongCortBigBigRx$CABIL_ID, response=LongCortBigBigRx$value, legend=T, col=c(1:5), xlab="Time", ylab="Log Cortisol", ylim=c(-1.5,0.5))
#+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))
quartz()
interaction.plot(LongCortBigRx$Time,trace.factor=LongCortBigRx$CABIL_ID, response=LongCortBigRx$value, legend=T, col=c(1:8), xlab="Time", ylab="Log Cortisol", ylim=c(-1.5,0.5))
#+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))
quartz()
interaction.plot(LongCortLittleRx$Time,trace.factor=LongCortLittleRx$CABIL_ID, response=LongCortLittleRx$value, legend=T, col="gray", xlab="Time", ylab="Log Cortisol", ylim=c(-1.5,0.5))
+axis(1, at=c(1,2), label=c("Before Mock Scan","After Mock Scan"))

# CortDxWideNoLog<-AllDataCortAnalysis[,c(1,2,8,9,77)]
# library(reshape)
# CortDxLongNoLog<-melt(CortDxWide, id=c("CABIL_ID","Dx", "CortLogDelta"))
# colnames(CortDxLongNoLog)[4]<-"Time"

# quartz()
# interaction.plot(CortDxLongNoLog$Time,trace.factor=CortDxLongNoLog$CABIL_ID, response=CortDxLongNoLog$value, legend=F, col=ifelse(abs(CortDxLongNoLog$CortLogDelta)>0.25, "red","black"))