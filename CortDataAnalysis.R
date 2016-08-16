setwd("/Users/abbiepopa/Documents/Lab/DPTB/Cortisol Analysis/cortdata from Elliot/Working")
Participants<-read.csv("PartList_DPTB.csv", na.strings=".")

RTTD<-read.csv("DPTB_RT_9-17-14_TDKoralyOut.csv",na.strings="NaN")
RTTD<-RTTD[,c(1,9,10)]

RT22q<-read.csv("DPTB_RT_9-17-14_22qKoralyOut.csv", na.strings="NaN")
RT22q<-RT22q[,c(1,9,10)]

RT<-rbind(RTTD,RT22q)


Cort<-read.csv("UniqueCort.csv")
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

AllDataCortAnalysisFull<-AllDataCortAnalysis


###note, use log10 

###Correlations, everyone together

CIx<-seq(-3,3, length.out=73)

#AllDataCortAnalysis<-na.omit(AllDataCortAnalysis)

	quartz()
	nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis$Spence.PARENT_PncAgr)),]
	plot(scale(nowData$LogPreCort), scale(as.numeric(nowData$Spence.PARENT_PncAgr)), xlab="LogPreCort", ylab="Spence.PARENT_PncAgr",col=ifelse(nowData$Dx=="22q","maroon2","blue"))
	abline(lm(scale(as.numeric(nowData$Spence.PARENT_PncAgr))~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(as.numeric(nowData$Spence.PARENT_PncAgr))~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(na.omit(nowData$LogPreCort)), decreasing=T),sort(na.omit(prd[,2]), decreasing=F),col="red", lty=2)
	  
	  lines(sort(scale(na.omit(nowData$LogPreCort)), decreasing=T),sort(na.omit(prd[,3]), decreasing=F),col="red", lty=2)
	  
	  #points(scale(nowData$LogPreCort),prd[,2],col="red")
	  
	 # points(scale(nowData$LogPreCort),prd[,3],col="red")
	 


	  


#76
for(i in 5:76){
	
	pretest<-cor.test(AllDataCortAnalysis$LogPreCort, AllDataCortAnalysis[,i])
	if(pretest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPreCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Dx=="22q","maroon2","blue"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)

	}
	posttest<-cor.test(AllDataCortAnalysis$LogPostCort, AllDataCortAnalysis[,i])
	if(posttest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$LogPostCort)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Dx=="22q","maroon2","blue"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
 
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	deltatest<-cor.test(AllDataCortAnalysis$CortDelta, AllDataCortAnalysis[,i])
	if(deltatest$p.value<0.05){
			nowData<-AllDataCortAnalysis[which(!is.na(AllDataCortAnalysis$CortDelta)&!is.na(AllDataCortAnalysis[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$Dx=="22q","maroon2","blue"))
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)

	}

	
	# logdeltatest<-cor.test(AllDataCortAnalysis$LogCortDelta, AllDataCortAnalysis[,i])
	# if(logdeltatest$p.value<0.05){
	# quartz()
	# plot(AllDataCortAnalysis$LogCortDelta, AllDataCortAnalysis[,i], main=paste("r=",round(logdeltatest$estimate, digits=2),"p=",round(logdeltatest$p.value, digits=2)), xlab="Log: Change in Cort (Post-Pre)", ylab=colnames(AllDataCortAnalysis)[i],col=ifelse(AllDataCortAnalysis$Dx=="22q","maroon2","blue"))
	# #abline(lm(AllDataCortAnalysis$LogCortDelta~AllDataCortAnalysis[,i]))
	# }
}

###change t-test, everone together
t.test(AllDataCortAnalysis$CortDelta)
#t.test(AllDataCortAnalysis$LogCortDelta)

###group t-tests
plot(AllDataCortAnalysis$Dx, AllDataCortAnalysis$LogPreCort)
t.test(AllDataCortAnalysis$LogPreCort[which(AllDataCortAnalysis$Dx=="22q")], AllDataCortAnalysis$LogPreCort[which(AllDataCortAnalysis$Dx=="TD")])
t.test(AllDataCortAnalysis$LogPostCort[which(AllDataCortAnalysis$Dx=="22q")], AllDataCortAnalysis$LogPostCort[which(AllDataCortAnalysis$Dx=="TD")])
t.test(AllDataCortAnalysis$CortDelta[which(AllDataCortAnalysis$Dx=="22q")], AllDataCortAnalysis$CortDelta[which(AllDataCortAnalysis$Dx=="TD")])
t.test(AllDataCortAnalysis$LogCortDelta[which(AllDataCortAnalysis$Dx=="22q")], AllDataCortAnalysis$LogCortDelta[which(AllDataCortAnalysis$Dx=="TD")])

###would change in looking times be a better measure than early and late looking times?

###Correlations, 22q
Cort22q<-subset(AllDataCortAnalysis, Dx=="22q")

for(i in 5:76){
	pretest<-cor.test(Cort22q$LogPreCort, Cort22q[,i])
	if(pretest$p.value<0.05){
		nowData<-Cort22q[which(!is.na(Cort22q$LogPreCort)&!is.na(Cort22q[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Dx=="22q","maroon2","blue"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	
	
	}
	posttest<-cor.test(Cort22q$LogPostCort, Cort22q[,i])
	if(posttest$p.value<0.05){
		nowData<-Cort22q[which(!is.na(Cort22q$LogPostCort)&!is.na(Cort22q[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Dx=="22q","maroon2","blue"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)

	}
	deltatest<-cor.test(Cort22q$CortDelta, Cort22q[,i])
	if(deltatest$p.value<0.05){
		nowData<-Cort22q[which(!is.na(Cort22q$CortDelta)&!is.na(Cort22q[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$Dx=="22q","maroon2","blue"))
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	# logdeltatest<-cor.test(Cort22q$LogCortDelta, Cort22q[,i])
	# if(logdeltatest$p.value<0.05){
	# quartz()
	# plot(Cort22q$LogCortDelta, Cort22q[,i], main=paste("r=",round(logdeltatest$estimate, digits=2),"p=",round(logdeltatest$p.value, digits=2)), xlab="Log: Change in Cort (Post-Pre)", ylab=colnames(Cort22q)[i],col=ifelse(Cort22q$Dx=="22q","maroon2","blue"))
	# #abline(lm(Cort22q$LogCortDelta~Cort22q[,i]))
	# }
}

###change t-test, 22q
t.test(Cort22q$CortDelta)
t.test(Cort22q$LogCortDelta)


###Correlations, TD
CortTD<-subset(AllDataCortAnalysis, Dx=="TD")

for(i in 5:76){
	pretest<-cor.test(CortTD$LogPreCort, CortTD[,i])
	if(pretest$p.value<0.05){
		nowData<-CortTD[which(!is.na(CortTD$LogPreCort)&!is.na(CortTD[,i])),]
	quartz()
	plot(scale(nowData$LogPreCort), scale(nowData[,i]), main=paste("r=",round(pretest$estimate, digits=2),"p=",round(pretest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPreCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Dx=="22q","maroon2","blue"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPreCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPreCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  
	  lines(sort(scale(nowData$LogPreCort), decreasing=(pretest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	
	
	}
	posttest<-cor.test(CortTD$LogPostCort, CortTD[,i])
	if(posttest$p.value<0.05){
		nowData<-CortTD[which(!is.na(CortTD$LogPostCort)&!is.na(CortTD[,i])),]
	quartz()
	plot(scale(nowData$LogPostCort), scale(nowData[,i]), main=paste("r=",round(posttest$estimate, digits=2),"p=",round(posttest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="LogPostCort", ylab=colnames(nowData)[i],col=ifelse(nowData$Dx=="22q","maroon2","blue"))
	abline(lm(scale(nowData[,i])~scale(nowData$LogPostCort)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$LogPostCort)), interval="confidence")
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  
	  lines(sort(scale(nowData$LogPostCort), decreasing=(posttest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)

	}
	deltatest<-cor.test(CortTD$CortDelta, CortTD[,i])
	if(deltatest$p.value<0.05){
		nowData<-CortTD[which(!is.na(CortTD$CortDelta)&!is.na(CortTD[,i])),]
	quartz()
	plot(scale(nowData$CortDelta), scale(nowData[,i]), main=paste("r=",round(deltatest$estimate, digits=2),"p=",round(deltatest$p.value, digits=2),"n=",dim(nowData)[[1]]), xlab="Change in Cort (Post-Pre)", ylab=colnames(nowData)[i],col=ifelse(nowData$Dx=="22q","maroon2","blue"))
	abline(lm(scale(nowData[,i])~scale(nowData$CortDelta)))
	  prd<-predict(lm(scale(nowData[,i])~scale(nowData$CortDelta)), interval="confidence")
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,2]),col="red", lty=2)
	  
	  lines(sort(scale(nowData$CortDelta), decreasing=(deltatest$estimate[[1]]<0)),sort(prd[,3]),col="red", lty=2)
	}
	# logdeltatest<-cor.test(Cort22q$LogCortDelta, Cort22q[,i])
	# if(logdeltatest$p.value<0.05){
	# quartz()
	# plot(Cort22q$LogCortDelta, Cort22q[,i], main=paste("r=",round(logdeltatest$estimate, digits=2),"p=",round(logdeltatest$p.value, digits=2)), xlab="Log: Change in Cort (Post-Pre)", ylab=colnames(Cort22q)[i],col=ifelse(Cort22q$Dx=="22q","maroon2","blue"))
	# #abline(lm(Cort22q$LogCortDelta~Cort22q[,i]))
	# }
}

###change t-test, TD
t.test(CortTD$CortDelta)
t.test(CortTD$LogCortDelta)


###Plots to understand cort values better
quartz()
ggplot(AllDataCortAnalysisFull, aes(factor(Dx),LogPreCort))+geom_violin(trim = F, aes(fill=factor(Dx)))+scale_fill_manual(values=c("maroon2","blue"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Diagnosis") + ylab("Log Pre Cortisol")
quartz()
ggplot(AllDataCortAnalysisFull, aes(factor(Dx),LogPostCort))+geom_violin(trim = F, aes(fill=factor(Dx)))+scale_fill_manual(values=c("maroon2","blue"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Diagnosis") + ylab("Log Post Cortisol")
quartz()
ggplot(AllDataCortAnalysisFull, aes(factor(Dx),CortDelta))+geom_violin(trim = F, aes(fill=factor(Dx)))+scale_fill_manual(values=c("maroon2","blue"))+geom_boxplot(width=0.1, fill="grey50") + xlab("Diagnosis") + ylab("Change in Cortisol")

CortDxWide<-AllDataCortAnalysis[,c(1,2,10,11,77)]
library(reshape)
CortDxLong<-melt(CortDxWide, id=c("CABIL_ID","Dx", "CortDelta"))
colnames(CortDxLong)[4]<-"Time"

# quartz()
# interaction.plot(CortDxLong$Time,trace.factor=CortDxLong$CABIL_ID, response=CortDxLong$value, legend=F, col=ifelse(CortDxLong$Dx=="22q","maroon2","blue"))

CortDxLongOrdered<-CortDxLong[order(CortDxLong$CABIL_ID),]
CortDxLongOrderedNoNA<-na.omit(CortDxLongOrdered)

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
	# MyRxColorVector[i]<-RxColoring(CortDxLongOrderedNoNA$CortDelta[i])
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

CortDeltaMean<-mean(CortDxLongOrderedNoNA$CortDelta)
CortDeltaSD<-sd(CortDxLongOrderedNoNA$CortDelta)
CortDelta1Up<-CortDeltaMean+CortDeltaSD
CortDelta1Down<-CortDeltaMean-CortDeltaSD
CortDelta2Up<-CortDeltaMean+(2*CortDeltaSD)
CortDelta2Down<-CortDeltaMean-(2*CortDeltaSD)

LongCortBigBigRx<-subset(CortDxLongOrderedNoNA, CortDelta>CortDelta2Up | CortDelta<CortDelta2Down)
LongCortBigRx<-subset(CortDxLongOrderedNoNA, CortDelta>CortDelta1Up & CortDelta<CortDelta2Up | CortDelta<CortDelta1Down & CortDelta>CortDelta2Down)
LongCortLittleRx<-subset(CortDxLongOrderedNoNA, CortDelta<CortDelta1Up & CortDelta>CortDelta1Down)

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
# CortDxLongNoLog<-melt(CortDxWide, id=c("CABIL_ID","Dx", "CortDelta"))
# colnames(CortDxLongNoLog)[4]<-"Time"

# quartz()
# interaction.plot(CortDxLongNoLog$Time,trace.factor=CortDxLongNoLog$CABIL_ID, response=CortDxLongNoLog$value, legend=F, col=ifelse(abs(CortDxLongNoLog$CortDelta)>0.25, "red","black"))