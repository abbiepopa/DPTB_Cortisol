#set working directory and important list of relevant participants
setwd("/Users/abbiepopa/Documents/Lab/DPTB/Cortisol Analysis/cortdata from Elliot/Working")
Participants<-read.csv("PartList_DPTB.csv", na.strings=".")

#import reaction time data for TD participants
RTTD<-read.csv("DPTB_RT_9-17-14_TDKoralyOut.csv",na.strings="NaN")
RTTD<-RTTD[,c(1,9,10)]

#import reaction time data for participants with 22q
RT22q<-read.csv("DPTB_RT_9-17-14_22qKoralyOut.csv", na.strings="NaN")
RT22q<-RT22q[,c(1,9,10)]

RT<-rbind(RTTD,RT22q)

#import cortisol data
Cort<-read.csv("CORT-Data-EAB_July2014.csv")
Cort<-Cort[,c(1, 3:6)]
colnames(Cort)[1]<-"CABIL_ID"

#import spence and ABAS data
SpenceABAS<-read.csv("Clustering_Database_8-8-14.csv", na.strings="-9999")

#import overall gaze data
EyeGazeOverall<-read.csv("AllData.csv", na.strings="NA")
EyeGazeOverall<-EyeGazeOverall[,2:11]
colnames(EyeGazeOverall)[1]<-"CABIL_ID"

#import time course gaze data
EyeGazeTimeCourse<-read.csv("AllDataTC25no676.csv", na.strings="NA")
EyeGazeTimeCourse<-EyeGazeTimeCourse[,2:24]
colnames(EyeGazeTimeCourse)[1]<-"CABIL_ID"

#import overall pupilometry for kids with 22q
PupilOverall22q<-read.csv("PupilPilot_6-5-14_22q.csv", na.strings=".")
PupilOverall22q<-PupilOverall22q[,c(1:7, 9:13)]
colnames(PupilOverall22q)[1]<-"CABIL_ID"

#import overall pupilometry to kids who are TD
PupilOverallTD<-read.csv("PupilPilot_6-5-14_TD.csv", na.strings=".")
PupilOverallTD<-PupilOverallTD[,c(1:7,9:13)]
colnames(PupilOverallTD)[1]<-"CABIL_ID"

#merge Dx groups
PupilOverall<-rbind(PupilOverall22q, PupilOverallTD)

#import time course pupilometry data for kids who are TD
PupilChangeTD<-read.csv("PupilChange_6-5-14_TD.csv",na.strings=".")
PupilChangeTD<-PupilChangeTD[,1:7]
colnames(PupilChangeTD)[1]<-"CABIL_ID"

#import time course pupilometry data for kids with 22q
PupilChange22q<-read.csv("PupilChange_6-5-14_22q.csv",na.strings=".")
PupilChange22q<-PupilChange22q[,1:7]
colnames(PupilChange22q)[1]<-"CABIL_ID"

#merge TC pupil data
PupilChange<-rbind(PupilChangeTD,PupilChange22q)

#merge data set
AllDataCortAnalysis<-merge(Participants, Cort, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, SpenceABAS, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, RT, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, EyeGazeOverall, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, EyeGazeTimeCourse, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, PupilOverall, all.x=T)
AllDataCortAnalysis<-merge(AllDataCortAnalysis, PupilChange, all.x=T)

#calculate change in cort, and change in log cort
AllDataCortAnalysis$CortDelta<-AllDataCortAnalysis$PostCORT-AllDataCortAnalysis$PreCORT
AllDataCortAnalysis$LogCortDelta<-log10(AllDataCortAnalysis$CortDelta+1)
AllDataCortAnalysis$CortLogDelta<-AllDataCortAnalysis$LogPostCort-AllDataCortAnalysis$LogPreCort

###clusters###
ClusAll<-read.csv("AllDataClus.csv")
Clus<-ClusAll[,c(2,12)]
colnames(Clus)[1]<-"CABIL_ID"

AllDataCortAnalysis<-merge(AllDataCortAnalysis, Clus, all.x=T)

### pos, neg, or no change in cort
BinAssigner<-function(x){
	if(is.na(x)){
		thebin<-NA
	}
	else if(x<(-0.10900)){
		thebin<-(-1)
	} else if (x>(0.10430)){
		thebin<-(1)
	}  else{
		thebin<-(0)
		}
		return(thebin)
}
for(i in 1:dim(AllDataCortAnalysis)[[1]]){
	AllDataCortAnalysis[i,c("CortDeltaBin")]<-BinAssigner(AllDataCortAnalysis[i,c("CortDelta")])
}

#Negative Only Graphs

NegChange<-subset(AllDataCortAnalysis, CortDeltaBin==-1)

for(i in 5:76){
	pretest<-cor.test(NegChange$LogPreCort, NegChange[,i])
	if(pretest$p.value<0.05){
			nowData<-NegChange[which(!is.na(NegChange$LogPreCort)&!is.na(NegChange[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col="darkblue",bg="blue", pch=25)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(NegChange$LogPostCort, NegChange[,i])
	if(posttest$p.value<0.05){
			nowData<-NegChange[which(!is.na(NegChange$LogPostCort)&!is.na(NegChange[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col="darkblue",bg="blue", pch=25)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(NegChange$CortDelta, NegChange[,i])
	if(deltatest$p.value<0.05){
			nowData<-NegChange[which(!is.na(NegChange$CortDelta)&!is.na(NegChange[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col="darkblue",bg="blue", pch=25)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#Negative Only Graphs

PosChange<-subset(AllDataCortAnalysis, CortDeltaBin==1)

for(i in 5:76){
	pretest<-cor.test(PosChange$LogPreCort, PosChange[,i])
	if(pretest$p.value<0.05){
			nowData<-PosChange[which(!is.na(PosChange$LogPreCort)&!is.na(PosChange[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col="darkgoldenrod4",bg="darkgoldenrod1", pch=24)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(PosChange$LogPostCort, PosChange[,i])
	if(posttest$p.value<0.05){
			nowData<-PosChange[which(!is.na(PosChange$LogPostCort)&!is.na(PosChange[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col="darkgoldenrod4",bg="darkgoldenrod1", pch=24)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(PosChange$CortDelta, PosChange[,i])
	if(deltatest$p.value<0.05){
			nowData<-PosChange[which(!is.na(PosChange$CortDelta)&!is.na(PosChange[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col="darkgoldenrod4",bg="darkgoldenrod1", pch=24)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#No Change Graphs

NoChange<-subset(AllDataCortAnalysis, CortDeltaBin==0)

for(i in 5:76){
	pretest<-cor.test(NoChange$LogPreCort, NoChange[,i])
	if(pretest$p.value<0.05){
			nowData<-NoChange[which(!is.na(NoChange$LogPreCort)&!is.na(NoChange[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col="darkgray",bg="lightgray", pch=22)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	posttest<-cor.test(NoChange$LogPostCort, NoChange[,i])
	if(posttest$p.value<0.05){
			nowData<-NoChange[which(!is.na(NoChange$LogPostCort)&!is.na(NoChange[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col="darkgray",bg="lightgray", pch=22)
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(NoChange$CortDelta, NoChange[,i])
	if(deltatest$p.value<0.05){
			nowData<-NoChange[which(!is.na(NoChange$CortDelta)&!is.na(NoChange[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col="darkgray",bg="lightgray", pch=22)
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
}

#Violins for Different Measures Based on ChangeInCort Bins
for (i in 5:76){
	thisfit<-Anova(lm(AllDataCortAnalysis[,i]~AllDataCortAnalysis$CortDeltaBin))
	if(thisfit[1,4]<0.05){
		quartz()
		print(ggplot(AllDataCortAnalysis, aes(factor(CortDeltaBin),AllDataCortAnalysis[,i]))+geom_violin(trim=F, aes(fill=factor(CortDeltaBin)))+scale_fill_manual(values=c("darkblue","lightgray","darkgoldenrod1"))+geom_boxplot(width=0.1,fill="grey50")+xlab("CortDeltaBin")+ylab(colnames(AllDataCortAnalysis)[i])+ggtitle(paste("p=",round(thisfit[1,4], digits=2))))
	}
}

###Follow-up Scatters to the Violins
library(stargazer)

#AngryEnd Angry

nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis$PercAngryEnd25_Angry)),]
fitblue<-lm(scale(nowData[which(nowData$CortDeltaBin<0),]$PercAngryEnd25_Angry)~scale(nowData[which(nowData$CortDeltaBin<0),]$CortDelta))
bluecor<-cor.test(nowData[which(nowData$CortDeltaBin<0),]$CortDelta, nowData[which(nowData$CortDeltaBin<0),]$PercAngryEnd25_Angry)
fitgray<-lm(scale(nowData[which(nowData$CortDeltaBin==0),]$PercAngryEnd25_Angry)~scale(nowData[which(nowData$CortDeltaBin==0),]$CortDelta))
graycor<-cor.test(nowData[which(nowData$CortDeltaBin==0),]$CortDelta, nowData[which(nowData$CortDeltaBin==0),]$PercAngryEnd25_Angry)
fitgold<-lm(scale(nowData[which(nowData$CortDeltaBin>0),]$PercAngryEnd25_Angry)~scale(nowData[which(nowData$CortDeltaBin>0),]$CortDelta))
goldcor<-cor.test(nowData[which(nowData$CortDeltaBin>0),]$CortDelta, nowData[which(nowData$CortDeltaBin>0),]$PercAngryEnd25_Angry)

quartz()
plot(scale(nowData[which(nowData$CortDeltaBin<0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBin<0),]$PercAngryEnd25_Angry), xlab="Change in Cort (Post-Pre)", ylab="Angry End - Angry", col="darkblue", bg="blue",
pch=25, main="CortBins")
abline(fitblue,col="blue")
points(scale(nowData[which(nowData$CortDeltaBin==0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBin==0),]$PercAngryEnd25_Angry), xlab="Change in Cort (Post-Pre)", ylab="Angry End - Angry", col="darkgray", bg="lightgray",
pch=22)
abline(fitgray, col="gray")
points(scale(nowData[which(nowData$CortDeltaBin>0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBi>0),]$PercAngryEnd25_Angry), xlab="Change in Cort (Post-Pre)", ylab="Angry End - Angry", col="darkgoldenrod4", bg="darkgoldenrod1",
pch=24)
abline(fitgold, col="darkgoldenrod3")

regtable<-matrix(nrow=3, ncol=3)
colnames(regtable)<-c("n","r","p")
row.names(regtable)<-c("Negative Change (blue)","NoChange (gray)","Positive Change (gold)")
regtable[1,1]<-dim(nowData[which(nowData$CortDeltaBin<0),])[[1]]
regtable[2,1]<-dim(nowData[which(nowData$CortDeltaBin==0),])[[1]]
regtable[3,1]<-dim(nowData[which(nowData$CortDeltaBin>0),])[[1]]

regtable[1,2]<-round(bluecor$estimate, digits=2)
regtable[2,2]<-round(graycor$estimate, digits=2)
regtable[3,2]<-round(goldcor$estimate, digits=2)

regtable[1,3]<-round(bluecor$p.value, digits=2)
regtable[2,3]<-round(graycor$p.value, digits=2)
regtable[3,3]<-round(goldcor$p.value, digits=2)

stargazer(regtable,title="Angry End - Angry",out="AngryEndAngry.html",summary=F)

#HappyChange Neutral

nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis$PercHappyChange25_Neutral)),]
fitblue<-lm(scale(nowData[which(nowData$CortDeltaBin<0),]$PercHappyChange25_Neutral)~scale(nowData[which(nowData$CortDeltaBin<0),]$CortDelta))
bluecor<-cor.test(nowData[which(nowData$CortDeltaBin<0),]$CortDelta, nowData[which(nowData$CortDeltaBin<0),]$PercHappyChange25_Neutral)
fitgray<-lm(scale(nowData[which(nowData$CortDeltaBin==0),]$PercHappyChange25_Neutral)~scale(nowData[which(nowData$CortDeltaBin==0),]$CortDelta))
graycor<-cor.test(nowData[which(nowData$CortDeltaBin==0),]$CortDelta, nowData[which(nowData$CortDeltaBin==0),]$PercHappyChange25_Neutral)
fitgold<-lm(scale(nowData[which(nowData$CortDeltaBin>0),]$PercHappyChange25_Neutral)~scale(nowData[which(nowData$CortDeltaBin>0),]$CortDelta))
goldcor<-cor.test(nowData[which(nowData$CortDeltaBin>0),]$CortDelta, nowData[which(nowData$CortDeltaBin>0),]$PercHappyChange25_Neutral)

quartz()
plot(scale(nowData[which(nowData$CortDeltaBin<0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBin<0),]$PercHappyChange25_Neutral), xlab="Change in Cort (Post-Pre)", ylab="Change in Neutral on Happy Trials", col="darkblue", bg="blue",
pch=25, main="CortBins")
abline(fitblue,col="blue")
points(scale(nowData[which(nowData$CortDeltaBin==0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBin==0),]$PercHappyChange25_Neutral), xlab="Change in Cort (Post-Pre)", ylab="Change in Neutral on Happy Trials", col="darkgray", bg="lightgray",
pch=22)
abline(fitgray, col="gray")
points(scale(nowData[which(nowData$CortDeltaBin>0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBi>0),]$PercHappyChange25_Neutral), xlab="Change in Cort (Post-Pre)", ylab="Change in Neutral on Happy Trials", col="darkgoldenrod4", bg="darkgoldenrod1",
pch=24)
abline(fitgold, col="darkgoldenrod3")

regtable<-matrix(nrow=3, ncol=3)
colnames(regtable)<-c("n","r","p")
row.names(regtable)<-c("Negative Change (blue)","NoChange (gray)","Positive Change (gold)")
regtable[1,1]<-dim(nowData[which(nowData$CortDeltaBin<0),])[[1]]
regtable[2,1]<-dim(nowData[which(nowData$CortDeltaBin==0),])[[1]]
regtable[3,1]<-dim(nowData[which(nowData$CortDeltaBin>0),])[[1]]

regtable[1,2]<-round(bluecor$estimate, digits=2)
regtable[2,2]<-round(graycor$estimate, digits=2)
regtable[3,2]<-round(goldcor$estimate, digits=2)

regtable[1,3]<-round(bluecor$p.value, digits=2)
regtable[2,3]<-round(graycor$p.value, digits=2)
regtable[3,3]<-round(goldcor$p.value, digits=2)

stargazer(regtable,title="Change in Neutral on Happy Trials",out="ChangeNeutralOnHappy.html",summary=F)

#HappyEnd Neutral

nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis$PercHappyEnd25_Neutral)),]
fitblue<-lm(scale(nowData[which(nowData$CortDeltaBin<0),]$PercHappyEnd25_Neutral)~scale(nowData[which(nowData$CortDeltaBin<0),]$CortDelta))
bluecor<-cor.test(nowData[which(nowData$CortDeltaBin<0),]$CortDelta, nowData[which(nowData$CortDeltaBin<0),]$PercHappyEnd25_Neutral)
fitgray<-lm(scale(nowData[which(nowData$CortDeltaBin==0),]$PercHappyEnd25_Neutral)~scale(nowData[which(nowData$CortDeltaBin==0),]$CortDelta))
graycor<-cor.test(nowData[which(nowData$CortDeltaBin==0),]$CortDelta, nowData[which(nowData$CortDeltaBin==0),]$PercHappyEnd25_Neutral)
fitgold<-lm(scale(nowData[which(nowData$CortDeltaBin>0),]$PercHappyEnd25_Neutral)~scale(nowData[which(nowData$CortDeltaBin>0),]$CortDelta))
goldcor<-cor.test(nowData[which(nowData$CortDeltaBin>0),]$CortDelta, nowData[which(nowData$CortDeltaBin>0),]$PercHappyEnd25_Neutral)

quartz()
plot(scale(nowData[which(nowData$CortDeltaBin<0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBin<0),]$PercHappyEnd25_Neutral), xlab="Change in Cort (Post-Pre)", ylab="Happy - Neutral End", col="darkblue", bg="blue",
pch=25, main="CortBins")
abline(fitblue,col="blue")
points(scale(nowData[which(nowData$CortDeltaBin==0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBin==0),]$PercHappyEnd25_Neutral), xlab="Change in Cort (Post-Pre)", ylab="Happy - Neutral End", col="darkgray", bg="lightgray",
pch=22)
abline(fitgray, col="gray")
points(scale(nowData[which(nowData$CortDeltaBin>0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBi>0),]$PercHappyEnd25_Neutral), xlab="Change in Cort (Post-Pre)", ylab="Happy - Neutral End", col="darkgoldenrod4", bg="darkgoldenrod1",
pch=24)
abline(fitgold, col="darkgoldenrod3")

regtable<-matrix(nrow=3, ncol=3)
colnames(regtable)<-c("n","r","p")
row.names(regtable)<-c("Negative Change (blue)","NoChange (gray)","Positive Change (gold)")
regtable[1,1]<-dim(nowData[which(nowData$CortDeltaBin<0),])[[1]]
regtable[2,1]<-dim(nowData[which(nowData$CortDeltaBin==0),])[[1]]
regtable[3,1]<-dim(nowData[which(nowData$CortDeltaBin>0),])[[1]]

regtable[1,2]<-round(bluecor$estimate, digits=2)
regtable[2,2]<-round(graycor$estimate, digits=2)
regtable[3,2]<-round(goldcor$estimate, digits=2)

regtable[1,3]<-round(bluecor$p.value, digits=2)
regtable[2,3]<-round(graycor$p.value, digits=2)
regtable[3,3]<-round(goldcor$p.value, digits=2)

stargazer(regtable,title="Happy - Neutral End",out="HappyNeutralEnd.html",summary=F)

#HappyEnd NonFace

nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis$PercHappyEnd25_NonFace)),]
fitblue<-lm(scale(nowData[which(nowData$CortDeltaBin<0),]$PercHappyEnd25_NonFace)~scale(nowData[which(nowData$CortDeltaBin<0),]$CortDelta))
bluecor<-cor.test(nowData[which(nowData$CortDeltaBin<0),]$CortDelta, nowData[which(nowData$CortDeltaBin<0),]$PercHappyEnd25_NonFace)
fitgray<-lm(scale(nowData[which(nowData$CortDeltaBin==0),]$PercHappyEnd25_NonFace)~scale(nowData[which(nowData$CortDeltaBin==0),]$CortDelta))
graycor<-cor.test(nowData[which(nowData$CortDeltaBin==0),]$CortDelta, nowData[which(nowData$CortDeltaBin==0),]$PercHappyEnd25_NonFace)
fitgold<-lm(scale(nowData[which(nowData$CortDeltaBin>0),]$PercHappyEnd25_NonFace)~scale(nowData[which(nowData$CortDeltaBin>0),]$CortDelta))
goldcor<-cor.test(nowData[which(nowData$CortDeltaBin>0),]$CortDelta, nowData[which(nowData$CortDeltaBin>0),]$PercHappyEnd25_NonFace)

quartz()
plot(scale(nowData[which(nowData$CortDeltaBin<0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBin<0),]$PercHappyEnd25_NonFace), xlab="Change in Cort (Post-Pre)", ylab="Happy - NonFace End", col="darkblue", bg="blue",
pch=25, main="CortBins")
abline(fitblue,col="blue")
points(scale(nowData[which(nowData$CortDeltaBin==0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBin==0),]$PercHappyEnd25_NonFace), xlab="Change in Cort (Post-Pre)", ylab="Happy - NonFace End", col="darkgray", bg="lightgray",
pch=22)
abline(fitgray, col="gray")
points(scale(nowData[which(nowData$CortDeltaBin>0),]$CortDelta), scale(nowData[which(nowData$CortDeltaBi>0),]$PercHappyEnd25_NonFace), xlab="Change in Cort (Post-Pre)", ylab="Happy - NonFace End", col="darkgoldenrod4", bg="darkgoldenrod1",
pch=24)
abline(fitgold, col="darkgoldenrod3")

regtable<-matrix(nrow=3, ncol=3)
colnames(regtable)<-c("n","r","p")
row.names(regtable)<-c("Negative Change (blue)","NoChange (gray)","Positive Change (gold)")
regtable[1,1]<-dim(nowData[which(nowData$CortDeltaBin<0),])[[1]]
regtable[2,1]<-dim(nowData[which(nowData$CortDeltaBin==0),])[[1]]
regtable[3,1]<-dim(nowData[which(nowData$CortDeltaBin>0),])[[1]]

regtable[1,2]<-round(bluecor$estimate, digits=2)
regtable[2,2]<-round(graycor$estimate, digits=2)
regtable[3,2]<-round(goldcor$estimate, digits=2)

regtable[1,3]<-round(bluecor$p.value, digits=2)
regtable[2,3]<-round(graycor$p.value, digits=2)
regtable[3,3]<-round(goldcor$p.value, digits=2)

stargazer(regtable,title="Happy - NonFace End",out="HappyNonFaceEnd.html",summary=F)

###tables about cortbins
for(i in c(2,80:84)){
	hereitis<-table(AllDataCortAnalysis[,i],AllDataCortAnalysis$CortDeltaBin)
	stargazer(as.data.frame.matrix(hereitis), title=colnames(AllDataCortAnalysis)[i],out=paste(colnames(AllDataCortAnalysis)[i],".html",sep=""), summary=F)
}
