setwd("/Users/abbiepopa/Documents/Lab/DPTB/Cortisol Analysis/cortdata from Elliot/Working")
Participants<-read.csv("PartList_DPTB.csv", na.strings=".")
library(ggplot2)

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

###Look at low General ABAS vs. high General ABAS
AllDataCortAnalysis$LogicABASGACImpaired<-(AllDataCortAnalysis$ABAS_GAC<76)

#Same Graphs
for(i in 5:76){
	pretest<-cor.test(AllDataCortAnalysis$LogPreCort, AllDataCortAnalysis[,i])
	if(pretest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$LogicABASGACImpaired, "red","forestgreen"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(AllDataCortAnalysis$LogPostCort, AllDataCortAnalysis[,i])
	if(posttest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPostCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$LogicABASGACImpaired, "red","forestgreen"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(AllDataCortAnalysis$CortDelta, AllDataCortAnalysis[,i])
	if(deltatest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$LogicABASGACImpaired, "red","forestgreen"))
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#Low Only
LowABASGAC<-subset(AllDataCortAnalysis, LogicABASGACImpaired)

for(i in 5:76){
	pretest<-cor.test(LowABASGAC$LogPreCort, LowABASGAC[,i])
	if(pretest$p.value<0.05){
			nowData<-LowABASGAC[which(!is.na(LowABASGAC$LogPreCort)&!is.na(LowABASGAC[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$LogicABASGACImpaired, "red","forestgreen"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(LowABASGAC$LogPostCort, LowABASGAC[,i])
	if(posttest$p.value<0.05){
			nowData<-LowABASGAC[which(!is.na(LowABASGAC$LogPostCort)&!is.na(LowABASGAC[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$LogicABASGACImpaired, "red","forestgreen"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(LowABASGAC$CortDelta, LowABASGAC[,i])
	if(deltatest$p.value<0.05){
			nowData<-LowABASGAC[which(!is.na(LowABASGAC$CortDelta)&!is.na(LowABASGAC[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$LogicABASGACImpaired, "red","forestgreen"))
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#High Only
HighABASGAC<-subset(AllDataCortAnalysis, !LogicABASGACImpaired)

for(i in 5:76){
	pretest<-cor.test(HighABASGAC$LogPreCort, HighABASGAC[,i])
	if(pretest$p.value<0.05){
			nowData<-HighABASGAC[which(!is.na(HighABASGAC$LogPreCort)&!is.na(HighABASGAC[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$LogicABASGACImpaired, "red","forestgreen"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(HighABASGAC$LogPostCort, HighABASGAC[,i])
	if(posttest$p.value<0.05){
			nowData<-HighABASGAC[which(!is.na(HighABASGAC$LogPostCort)&!is.na(HighABASGAC[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$LogicABASGACImpaired, "red","forestgreen"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(HighABASGAC$CortDelta, HighABASGAC[,i])
	if(deltatest$p.value<0.05){
			nowData<-HighABASGAC[which(!is.na(HighABASGAC$CortDelta)&!is.na(HighABASGAC[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$LogicABASGACImpaired, "red","forestgreen"))
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

###Look at low Total Spence vs high Total Spence
AllDataCortAnalysis$Spnc_Eleveated<-(AllDataCortAnalysis$Spnc_Ttl_Parent>59)

#Same Graphs
for(i in 5:76){
	pretest<-cor.test(AllDataCortAnalysis$LogPreCort, AllDataCortAnalysis[,i])
	if(pretest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Spnc_Eleveated, "red","forestgreen"), pch=8)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(AllDataCortAnalysis$LogPostCort, AllDataCortAnalysis[,i])
	if(posttest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPostCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Spnc_Eleveated, "red","forestgreen"), pch=8)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(AllDataCortAnalysis$CortDelta, AllDataCortAnalysis[,i])
	if(deltatest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$Spnc_Eleveated, "red","forestgreen"), pch=8)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#High Spence
HighSpncTtl<-subset(AllDataCortAnalysis, Spnc_Eleveated)

for(i in 5:76){
	pretest<-cor.test(HighSpncTtl$LogPreCort, HighSpncTtl[,i])
	if(pretest$p.value<0.05){
			nowData<-HighSpncTtl[which(!is.na(HighSpncTtl$LogPreCort)&!is.na(HighSpncTtl[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Spnc_Eleveated, "red","forestgreen"), pch=8)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(HighSpncTtl$LogPostCort, HighSpncTtl[,i])
	if(posttest$p.value<0.05){
			nowData<-HighSpncTtl[which(!is.na(HighSpncTtl$LogPostCort)&!is.na(HighSpncTtl[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Spnc_Eleveated, "red","forestgreen"), pch=8)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(HighSpncTtl$CortDelta, HighSpncTtl[,i])
	if(deltatest$p.value<0.05){
			nowData<-HighSpncTtl[which(!is.na(HighSpncTtl$CortDelta)&!is.na(HighSpncTtl[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$Spnc_Eleveated, "red","forestgreen"), pch=8)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#Low Spence
LowSpncTtl<-subset(AllDataCortAnalysis, !Spnc_Eleveated)

for(i in 5:76){
	pretest<-cor.test(LowSpncTtl$LogPreCort, LowSpncTtl[,i])
	if(pretest$p.value<0.05){
			nowData<-LowSpncTtl[which(!is.na(LowSpncTtl$LogPreCort)&!is.na(LowSpncTtl[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Spnc_Eleveated, "red","forestgreen"), pch=8)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(LowSpncTtl$LogPostCort, LowSpncTtl[,i])
	if(posttest$p.value<0.05){
			nowData<-LowSpncTtl[which(!is.na(LowSpncTtl$LogPostCort)&!is.na(LowSpncTtl[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Spnc_Eleveated, "red","forestgreen"), pch=8)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(LowSpncTtl$CortDelta, LowSpncTtl[,i])
	if(deltatest$p.value<0.05){
			nowData<-LowSpncTtl[which(!is.na(LowSpncTtl$CortDelta)&!is.na(LowSpncTtl[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$Spnc_Eleveated, "red","forestgreen"), pch=8)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

###Look at low physical injury vs high physical injury
AllDataCortAnalysis$PhysInjElevated<-(AllDataCortAnalysis$Spence.PARENT_PhysInj>59)

#Same Graphs
for(i in 5:76){
	pretest<-cor.test(AllDataCortAnalysis$LogPreCort, AllDataCortAnalysis[,i])
	if(pretest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$PhysInjElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(AllDataCortAnalysis$LogPostCort, AllDataCortAnalysis[,i])
	if(posttest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPostCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$PhysInjElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(AllDataCortAnalysis$CortDelta, AllDataCortAnalysis[,i])
	if(deltatest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$PhysInjElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#High Phys Inj
HighPhysInj<-subset(AllDataCortAnalysis, PhysInjElevated)

for(i in 5:76){
	pretest<-cor.test(HighPhysInj$LogPreCort, HighPhysInj[,i])
	if(pretest$p.value<0.05){
			nowData<-HighPhysInj[which(!is.na(HighPhysInj$LogPreCort)&!is.na(HighPhysInj[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$PhysInjElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(HighPhysInj$LogPostCort, HighPhysInj[,i])
	if(posttest$p.value<0.05){
			nowData<-HighPhysInj[which(!is.na(HighPhysInj$LogPostCort)&!is.na(HighPhysInj[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$PhysInjElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(HighPhysInj$CortDelta, HighPhysInj[,i])
	if(deltatest$p.value<0.05){
			nowData<-HighPhysInj[which(!is.na(HighPhysInj$CortDelta)&!is.na(HighPhysInj[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$PhysInjElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#Lwo Phys Inj
LowPhysInj<-subset(AllDataCortAnalysis, !PhysInjElevated)

for(i in 5:76){
	pretest<-cor.test(LowPhysInj$LogPreCort, LowPhysInj[,i])
	if(pretest$p.value<0.05){
			nowData<-LowPhysInj[which(!is.na(LowPhysInj$LogPreCort)&!is.na(LowPhysInj[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$PhysInjElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(LowPhysInj$LogPostCort, LowPhysInj[,i])
	if(posttest$p.value<0.05){
			nowData<-LowPhysInj[which(!is.na(LowPhysInj$LogPostCort)&!is.na(LowPhysInj[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$PhysInjElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(LowPhysInj$CortDelta, LowPhysInj[,i])
	if(deltatest$p.value<0.05){
			nowData<-LowPhysInj[which(!is.na(LowPhysInj$CortDelta)&!is.na(LowPhysInj[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$PhysInjElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

###Look at low separation anxiety vs high separation anxiety
AllDataCortAnalysis$SepAnxElevated<-(AllDataCortAnalysis$Spence.PARENT_SepAnx>59)

#Same Graphs
for(i in 5:76){
	pretest<-cor.test(AllDataCortAnalysis$LogPreCort, AllDataCortAnalysis[,i])
	if(pretest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$SepAnxElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(AllDataCortAnalysis$LogPostCort, AllDataCortAnalysis[,i])
	if(posttest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPostCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$SepAnxElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(AllDataCortAnalysis$CortDelta, AllDataCortAnalysis[,i])
	if(deltatest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$SepAnxElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#High Sep Anx
HighSepAnx<-subset(AllDataCortAnalysis, SepAnxElevated)

for(i in 5:76){
	pretest<-cor.test(HighSepAnx$LogPreCort, HighSepAnx[,i])
	if(pretest$p.value<0.05){
			nowData<-HighSepAnx[which(!is.na(HighSepAnx$LogPreCort)&!is.na(HighSepAnx[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$SepAnxElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(HighSepAnx$LogPostCort, HighSepAnx[,i])
	if(posttest$p.value<0.05){
			nowData<-HighSepAnx[which(!is.na(HighSepAnx$LogPostCort)&!is.na(HighSepAnx[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$SepAnxElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(HighSepAnx$CortDelta, HighSepAnx[,i])
	if(deltatest$p.value<0.05){
			nowData<-HighSepAnx[which(!is.na(HighSepAnx$CortDelta)&!is.na(HighSepAnx[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$SepAnxElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#Low Sep Anx
LowSepAnx<-subset(AllDataCortAnalysis, !SepAnxElevated)

for(i in 5:76){
	pretest<-cor.test(LowSepAnx$LogPreCort, LowSepAnx[,i])
	if(pretest$p.value<0.05){
			nowData<-LowSepAnx[which(!is.na(LowSepAnx$LogPreCort)&!is.na(LowSepAnx[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$SepAnxElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(LowSepAnx$LogPostCort, LowSepAnx[,i])
	if(posttest$p.value<0.05){
			nowData<-LowSepAnx[which(!is.na(LowSepAnx$LogPostCort)&!is.na(LowSepAnx[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$SepAnxElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(LowSepAnx$CortDelta, LowSepAnx[,i])
	if(deltatest$p.value<0.05){
			nowData<-LowSepAnx[which(!is.na(LowSepAnx$CortDelta)&!is.na(LowSepAnx[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$SepAnxElevated, "red","forestgreen"), pch=2)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

###Look at all SCAS correlations for only those who are impaired on that measure

for(i in c(5,6,12:17)){
		preData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]>59),]
		pretest<-cor.test(preData$LogPreCort, preData[,i])
	if(pretest$p.value<0.05){
	quartz()
	plot(scale(preData$LogPreCort), scale(preData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(preData)[[1]]), xlab="LogPreCort", ylab=colnames(preData)[i], col="red", pch=12)
	abline(lm(scale(preData[,i])~scale(preData$LogPreCort)))
	  prd<-predict(lm(scale(preData[,i])~scale(preData$LogPreCort)), interval="confidence")
	  lines(sort(scale(preData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(preData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	postData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPostCort)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]>59),]
	posttest<-cor.test(postData$LogPostCort, postData[,i])
	if(posttest$p.value<0.05){
			
	quartz()
	plot(scale(postData$LogPostCort), scale(postData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(postData)[[1]]), xlab="LogPostCort", ylab=colnames(postData)[i],col="red", pch=12)
	abline(lm(scale(postData[,i])~scale(postData$LogPostCort)))
	  prd<-predict(lm(scale(postData[,i])~scale(postData$LogPostCort)), interval="confidence")
	  lines(sort(scale(postData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(postData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltaData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]>59),]
	deltatest<-cor.test(deltaData$CortDelta, deltaData[,i])
	if(deltatest$p.value<0.05){
			
	quartz()
	plot(scale(deltaData$CortDelta), scale(deltaData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(deltaData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(deltaData)[i],col="red", pch=12)
	abline(lm(scale(deltaData[,i])~scale(deltaData$CortDelta)))
	  prd<-predict(lm(scale(deltaData[,i])~scale(deltaData$CortDelta)), interval="confidence")
	  lines(sort(scale(deltaData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(deltaData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	}

###Unelevated Spence
for(i in c(5,6,12:17)){
		preData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]<60),]
		pretest<-cor.test(preData$LogPreCort, preData[,i])
	if(pretest$p.value<0.05){
	quartz()
	plot(scale(preData$LogPreCort), scale(preData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(preData)[[1]]), xlab="LogPreCort", ylab=colnames(preData)[i], col="forestgreen", pch=12)
	abline(lm(scale(preData[,i])~scale(preData$LogPreCort)))
	  prd<-predict(lm(scale(preData[,i])~scale(preData$LogPreCort)), interval="confidence")
	  lines(sort(scale(preData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(preData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	postData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPostCort)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]<60),]
	posttest<-cor.test(postData$LogPostCort, postData[,i])
	if(posttest$p.value<0.05){
			
	quartz()
	plot(scale(postData$LogPostCort), scale(postData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(postData)[[1]]), xlab="LogPostCort", ylab=colnames(postData)[i],col="forestgreen", pch=12)
	abline(lm(scale(postData[,i])~scale(postData$LogPostCort)))
	  prd<-predict(lm(scale(postData[,i])~scale(postData$LogPostCort)), interval="confidence")
	  lines(sort(scale(postData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(postData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltaData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]<60),]
	deltatest<-cor.test(deltaData$CortDelta, deltaData[,i])
	if(deltatest$p.value<0.05){
			
	quartz()
	plot(scale(deltaData$CortDelta), scale(deltaData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(deltaData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(deltaData)[i],col="forestgreen", pch=12)
	abline(lm(scale(deltaData[,i])~scale(deltaData$CortDelta)))
	  prd<-predict(lm(scale(deltaData[,i])~scale(deltaData$CortDelta)), interval="confidence")
	  lines(sort(scale(deltaData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(deltaData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

###Look at all the ABAS correlations for only those who are impaired on that measure

for(i in c(18:26)){
		preData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]<8),]
		pretest<-cor.test(preData$LogPreCort, preData[,i])
	if(pretest$p.value<0.05){
	quartz()
	plot(scale(preData$LogPreCort), scale(preData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(preData)[[1]]), xlab="LogPreCort", ylab=colnames(preData)[i], col="red", pch=9)
	abline(lm(scale(preData[,i])~scale(preData$LogPreCort)))
	  prd<-predict(lm(scale(preData[,i])~scale(preData$LogPreCort)), interval="confidence")
	  lines(sort(scale(preData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(preData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	postData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPostCort)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]<8),]
	posttest<-cor.test(postData$LogPostCort, postData[,i])
	if(posttest$p.value<0.05){
			
	quartz()
	plot(scale(postData$LogPostCort), scale(postData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(postData)[[1]]), xlab="LogPostCort", ylab=colnames(postData)[i],col="red", pch=9)
	abline(lm(scale(postData[,i])~scale(postData$LogPostCort)))
	  prd<-predict(lm(scale(postData[,i])~scale(postData$LogPostCort)), interval="confidence")
	  lines(sort(scale(postData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(postData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltaData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]<8),]
	deltatest<-cor.test(deltaData$CortDelta, deltaData[,i])
	if(deltatest$p.value<0.05){
			
	quartz()
	plot(scale(deltaData$CortDelta), scale(deltaData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(deltaData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(deltaData)[i],col="red", pch=9)
	abline(lm(scale(deltaData[,i])~scale(deltaData$CortDelta)))
	  prd<-predict(lm(scale(deltaData[,i])~scale(deltaData$CortDelta)), interval="confidence")
	  lines(sort(scale(deltaData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(deltaData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

for(i in c(18:26)){
		preData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]>7),]
		pretest<-cor.test(preData$LogPreCort, preData[,i])
	if(pretest$p.value<0.05){
	quartz()
	plot(scale(preData$LogPreCort), scale(preData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(preData)[[1]]), xlab="LogPreCort", ylab=colnames(preData)[i], col="forestgreen", pch=9)
	abline(lm(scale(preData[,i])~scale(preData$LogPreCort)))
	  prd<-predict(lm(scale(preData[,i])~scale(preData$LogPreCort)), interval="confidence")
	  lines(sort(scale(preData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(preData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	postData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPostCort)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]>7),]
	posttest<-cor.test(postData$LogPostCort, postData[,i])
	if(posttest$p.value<0.05){
			
	quartz()
	plot(scale(postData$LogPostCort), scale(postData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(postData)[[1]]), xlab="LogPostCort", ylab=colnames(postData)[i],col="forestgreen", pch=9)
	abline(lm(scale(postData[,i])~scale(postData$LogPostCort)))
	  prd<-predict(lm(scale(postData[,i])~scale(postData$LogPostCort)), interval="confidence")
	  lines(sort(scale(postData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(postData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltaData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis[,i])&AllDataCortAnalysis[,i]>7),]
	deltatest<-cor.test(deltaData$CortDelta, deltaData[,i])
	if(deltatest$p.value<0.05){
			
	quartz()
	plot(scale(deltaData$CortDelta), scale(deltaData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(deltaData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(deltaData)[i],col="forestgreen", pch=9)
	abline(lm(scale(deltaData[,i])~scale(deltaData$CortDelta)))
	  prd<-predict(lm(scale(deltaData[,i])~scale(deltaData$CortDelta)), interval="confidence")
	  lines(sort(scale(deltaData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(deltaData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

###t.tests & violins

#ABAS GAC 
t.test(AllDataCortAnalysis$LogPreCort[which(AllDataCortAnalysis$LogicABASGACImpaired)], AllDataCortAnalysis$LogPreCort[which(!AllDataCortAnalysis$LogicABASGACImpaired)])
t.test(AllDataCortAnalysis$LogPostCort[which(AllDataCortAnalysis$LogicABASGACImpaired)], AllDataCortAnalysis$LogPostCort[which(!AllDataCortAnalysis$LogicABASGACImpaired)])
t.test(AllDataCortAnalysis$CortDelta[which(AllDataCortAnalysis$LogicABASGACImpaired)], AllDataCortAnalysis$CortDelta[which(!AllDataCortAnalysis$LogicABASGACImpaired)])


quartz()
ggplot(AllDataCortAnalysis, aes(factor(LogicABASGACImpaired),LogPreCort))+geom_violin(trim = F, aes(fill=factor(LogicABASGACImpaired)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("ABAS GAC impaired?") + ylab("Log Pre Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(LogicABASGACImpaired),LogPostCort))+geom_violin(trim = F, aes(fill=factor(LogicABASGACImpaired)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("ABAS GAC impaired?") + ylab("Log Post Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(LogicABASGACImpaired),CortDelta))+geom_violin(trim = F, aes(fill=factor(LogicABASGACImpaired)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("ABAS GAC impaired?") + ylab("Change in Cortisol")

#Spence Parent Total
t.test(AllDataCortAnalysis$LogPreCort[which(AllDataCortAnalysis$Spnc_Eleveated)], AllDataCortAnalysis$LogPreCort[which(!AllDataCortAnalysis$Spnc_Eleveated)])
t.test(AllDataCortAnalysis$LogPostCort[which(AllDataCortAnalysis$Spnc_Eleveated)], AllDataCortAnalysis$LogPostCort[which(!AllDataCortAnalysis$Spnc_Eleveated)])
t.test(AllDataCortAnalysis$CortDelta[which(AllDataCortAnalysis$Spnc_Eleveated)], AllDataCortAnalysis$CortDelta[which(!AllDataCortAnalysis$Spnc_Eleveated)])

quartz()
ggplot(AllDataCortAnalysis, aes(factor(Spnc_Eleveated),LogPreCort))+geom_violin(trim = F, aes(fill=factor(Spnc_Eleveated)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Spence Parent Total Elevated?") + ylab("Log Pre Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(Spnc_Eleveated),LogPostCort))+geom_violin(trim = F, aes(fill=factor(Spnc_Eleveated)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Spence Parent Total Elevated?") + ylab("Log Post Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(Spnc_Eleveated),CortDelta))+geom_violin(trim = F, aes(fill=factor(Spnc_Eleveated)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Spence Parent Total Elevated?") + ylab("Change in Cortisol")

#Spence Phys Inj
t.test(AllDataCortAnalysis$LogPreCort[which(AllDataCortAnalysis$PhysInjElevated)], AllDataCortAnalysis$LogPreCort[which(!AllDataCortAnalysis$PhysInjElevated)])
t.test(AllDataCortAnalysis$LogPostCort[which(AllDataCortAnalysis$PhysInjElevated)], AllDataCortAnalysis$LogPostCort[which(!AllDataCortAnalysis$PhysInjElevated)])
t.test(AllDataCortAnalysis$CortDelta[which(AllDataCortAnalysis$PhysInjElevated)], AllDataCortAnalysis$CortDelta[which(!AllDataCortAnalysis$PhysInjElevated)])

quartz()
ggplot(AllDataCortAnalysis, aes(factor(PhysInjElevated),LogPreCort))+geom_violin(trim = F, aes(fill=factor(PhysInjElevated)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Spence Parent Physical Injury Elevated?") + ylab("Log Pre Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(PhysInjElevated),LogPostCort))+geom_violin(trim = F, aes(fill=factor(PhysInjElevated)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Spence Parent Physical Injury Elevated?") + ylab("Log Post Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(PhysInjElevated),CortDelta))+geom_violin(trim = F, aes(fill=factor(PhysInjElevated)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Spence Parent Physical Injury Elevated?") + ylab("Change in Cortisol")

#Spence Sep Anx
t.test(AllDataCortAnalysis$LogPreCort[which(AllDataCortAnalysis$SepAnxElevated)], AllDataCortAnalysis$LogPreCort[which(!AllDataCortAnalysis$SepAnxElevated)])
t.test(AllDataCortAnalysis$LogPostCort[which(AllDataCortAnalysis$SepAnxElevated)], AllDataCortAnalysis$LogPostCort[which(!AllDataCortAnalysis$SepAnxElevated)])
t.test(AllDataCortAnalysis$CortDelta[which(AllDataCortAnalysis$SepAnxElevated)], AllDataCortAnalysis$CortDelta[which(!AllDataCortAnalysis$SepAnxElevated)])

quartz()
ggplot(AllDataCortAnalysis, aes(factor(SepAnxElevated),LogPreCort))+geom_violin(trim = F, aes(fill=factor(SepAnxElevated)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Spence Parent Separation Anxiety Elevated?") + ylab("Log Pre Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(SepAnxElevated),LogPostCort))+geom_violin(trim = F, aes(fill=factor(SepAnxElevated)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Spence Parent Separation Anxiety Elevated?") + ylab("Log Post Cortisol")
quartz()
ggplot(AllDataCortAnalysis, aes(factor(SepAnxElevated),CortDelta))+geom_violin(trim = F, aes(fill=factor(SepAnxElevated)))+scale_fill_manual(values=c("forestgreen","red"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Spence Parent Separation Anxiety Elevated?") + ylab("Change in Cortisol")


