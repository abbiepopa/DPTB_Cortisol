opar<-par()


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

### pow, neg, or no change in cort
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

PosChange<-subset(AllDataCortAnalysis, CortDeltaBin==1)

NoChange<-subset(AllDataCortAnalysis, CortDeltaBin==0)

fit1<-lm(NegChange$PercAngryBegin25_Angry~NegChange$LogPreCort)
prd1<-predict(fit1,interval="confidence")

par(cex.axis=1.8)
par(cex.main=2)
par(cex.lab=1.8)
par(mar=c(5.1,5,4.6,2.1))
plot(NegChange[which(NegChange$Dx =="22q"),"LogPreCort"], NegChange[which(NegChange$Dx =="22q"),"PercAngryBegin25_Angry"],col="Dark Blue", pch=20, cex=2, ylim=c(0,.5), main = "Angry Face Viewing vs. Cortisol \nin Reducers", xlab="Log Pre Cortisol", ylab="\nRatio Angry Face Viewing - Early Trials")
points(NegChange[which(NegChange$Dx =="TD"),"LogPreCort"], NegChange[which(NegChange$Dx =="TD"),"PercAngryBegin25_Angry"], col="DarkBlue", pch=17, cex=2)
abline(fit1, lwd=2)
lines(sort((NegChange[which(!is.na(NegChange$PercAngryBegin25_Angry)),"LogPreCort"])), sort(na.omit(prd1[,2])), col="red", lty=2, lwd=2)
lines(sort((NegChange[which(!is.na(NegChange$PercAngryBegin25_Angry)),"LogPreCort"])), sort(na.omit(prd1[,3])), col="red", lty=2, lwd=2)

plot(c(1,10),c(1,10))
legend(2,8,c("Reducers with 22q", "TD Reducers", "Fit Line r=0.61", "p=0.04", "95% CI"), col=c("dark blue", "dark blue", "black",NA, "red"), pch=c(20,17, NA, NA,NA), lty=c(NA, NA, 1, NA, 2), lwd=c(NA, NA, 2,NA, 2), cex=2)



fit2<-lm(NegChange$PercAngryBegin25_NonFace~NegChange$LogPreCort)
prd2<-predict(fit2,interval="confidence")

plot(NegChange[which(NegChange$Dx =="22q"),"LogPreCort"], NegChange[which(NegChange$Dx =="22q"),"PercAngryBegin25_NonFace"],col="Dark Blue", pch=20, cex=2, ylim=c(0.1,.9), main = "NonFace Viewing vs. Cortisol \nin Reducers", xlab="Log Pre Cortisol", ylab="Ratio NonFace Viewing - Early Trials")
points(NegChange[which(NegChange$Dx =="TD"),"LogPreCort"], NegChange[which(NegChange$Dx =="TD"),"PercAngryBegin25_NonFace"], col="DarkBlue", pch=17, cex=2)
abline(fit2, lwd=2)
lines(sort((NegChange[which(!is.na(NegChange$PercAngryBegin25_NonFace)),"LogPreCort"])), sort(na.omit(prd2[,2]), decreasing=T), col="red", lty=2, lwd=2)
lines(sort((NegChange[which(!is.na(NegChange$PercAngryBegin25_NonFace)),"LogPreCort"])), sort(na.omit(prd2[,3]), decreasing=T), col="red", lty=2, lwd=2)

plot(c(1,10),c(1,10))
legend(2,8,c("Reducers with 22q", "TD Reducers", "Fit Line r=-0.78", "p=0.003", "95% CI"), col=c("dark blue", "dark blue", "black",NA, "red"), pch=c(20,17, NA, NA,NA), lty=c(NA, NA, 1, NA, 2), lwd=c(NA, NA, 2,NA, 2), cex=2)

fit1_inc<-lm(PosChange$PercAngryBegin25_Angry~PosChange$LogPreCort)
fit1_main<-lm(NoChange$PercAngryBegin25_Angry~NoChange$LogPreCort)

plot(NegChange[which(NegChange$Dx =="22q"),"LogPreCort"], NegChange[which(NegChange$Dx =="22q"),"PercAngryBegin25_Angry"],col="Dark Blue", pch=20, cex=2, ylim=c(0,0.7), main = "Angry Face Viewing vs. Cortisol \nin All Participants", xlab="Log Pre Cortisol", ylab="Ratio Angry Face Viewing - Early Trials")
points(NegChange[which(NegChange$Dx =="TD"),"LogPreCort"], NegChange[which(NegChange$Dx =="TD"),"PercAngryBegin25_Angry"], col="DarkBlue", pch=17, cex=2)
abline(fit1, lwd=2, col="Dark Blue")


points(PosChange[which(PosChange$Dx =="22q"),"LogPreCort"], PosChange[which(PosChange$Dx =="22q"),"PercAngryBegin25_Angry"],col="Dark Red", pch=20, cex=2)
points(PosChange[which(PosChange$Dx =="TD"),"LogPreCort"], PosChange[which(PosChange$Dx =="TD"),"PercAngryBegin25_Angry"], col="Dark Red", pch=17, cex=2)
abline(fit1_inc, lwd=2, col="Dark Red")

points(NoChange[which(NoChange$Dx =="22q"),"LogPreCort"], NoChange[which(NoChange$Dx =="22q"),"PercAngryBegin25_Angry"],col="darkgoldenrod1", pch=20, cex=2)
points(NoChange[which(NoChange$Dx =="TD"),"LogPreCort"], NoChange[which(NoChange$Dx =="TD"),"PercAngryBegin25_Angry"], col="darkgoldenrod1", pch=17, cex=2)
abline(fit1_main, lwd=2, col="darkgoldenrod1")

fit2_inc<-lm(PosChange$PercAngryBegin25_NonFace~PosChange$LogPreCort)
fit2_main<-lm(NoChange$PercAngryBegin25_NonFace~NoChange$LogPreCort)

plot(NegChange[which(NegChange$Dx =="22q"),"LogPreCort"], NegChange[which(NegChange$Dx =="22q"),"PercAngryBegin25_NonFace"],col="Dark Blue", pch=20, cex=2, ylim=c(0,1), main = "NonFace Viewing vs. Cortisol \nin All Participants", xlab="Log Pre Cortisol", ylab="Ratio NonFace Viewing - Early Trials")
points(NegChange[which(NegChange$Dx =="TD"),"LogPreCort"], NegChange[which(NegChange$Dx =="TD"),"PercAngryBegin25_NonFace"], col="DarkBlue", pch=17, cex=2)
abline(fit2, lwd=2, col="Dark Blue")


points(PosChange[which(PosChange$Dx =="22q"),"LogPreCort"], PosChange[which(PosChange$Dx =="22q"),"PercAngryBegin25_NonFace"],col="Dark Red", pch=20, cex=2)
points(PosChange[which(PosChange$Dx =="TD"),"LogPreCort"], PosChange[which(PosChange$Dx =="TD"),"PercAngryBegin25_NonFace"], col="Dark Red", pch=17, cex=2)
abline(fit2_inc, lwd=2, col="Dark Red")

points(NoChange[which(NoChange$Dx =="22q"),"LogPreCort"], NoChange[which(NoChange$Dx =="22q"),"PercAngryBegin25_NonFace"],col="darkgoldenrod1", pch=20, cex=2)
points(NoChange[which(NoChange$Dx =="TD"),"LogPreCort"], NoChange[which(NoChange$Dx =="TD"),"PercAngryBegin25_NonFace"], col="darkgoldenrod1", pch=17, cex=2)
abline(fit2_main, lwd=2, col="darkgoldenrod1")

plot(c(1,10),c(1,10))
legend(2,9,c("22q","TD"),pch=c(20,17),cex=2)

###Legend for ANGRY all
plot(c(1,10),c(1,10))
legend(1,9,c("Reducers: r=0.61 p=0.04","Maintainers: r=0.21 p=0.41","Increasers: r=-0.27 p=0.44"), col=c("Dark Blue","darkgoldenrod1","Dark Red"), lwd=c(2,2,2), cex=2)

###Legend for NONFACE all
plot(c(0.5,10),c(1,10))
legend(0.3,9,c("Reducers: r=-0.78 p=0.003", "Maintainers: r=-0.24 p=0.35", "Increasers: r=0.38 p=0.29"), col=c("Dark Blue","darkgoldenrod1","Dark Red"), lwd=c(2,2,2), cex=2)

###Early vs. Late
fit1_late<-lm(NegChange$PercAngryEnd25_Angry~NegChange$LogPostCort)

plot(NegChange[which(NegChange$Dx =="22q"),"LogPreCort"], NegChange[which(NegChange$Dx =="22q"),"PercAngryBegin25_Angry"],col="Dark Blue", pch=20, cex=2, ylim=c(0,0.6), xlim=c(-.5,.4), main = "Early and Late Gaze in Reducers  \nAngry Faces", xlab="Log Cortisol", ylab="Ratio Angry Face Viewing")
points(NegChange[which(NegChange$Dx =="TD"),"LogPreCort"], NegChange[which(NegChange$Dx =="TD"),"PercAngryBegin25_Angry"], col="DarkBlue", pch=17, cex=2)
abline(fit1, lwd=2, col="Dark Blue")

points(NegChange[which(NegChange$Dx =="22q"),"LogPreCort"], NegChange[which(NegChange$Dx =="22q"),"PercAngryEnd25_Angry"], col="Light Blue", pch=17, cex=2)
points(NegChange[which(NegChange$Dx =="TD"),"LogPreCort"], NegChange[which(NegChange$Dx =="TD"),"PercAngryEnd25_Angry"], col="Light Blue", pch=17, cex=2)
abline(fit1_late, lwd=2, col="Light Blue")


plot(c(1,10),c(1,10))
legend(2, 9, c("Early", "r=0.61", "p=0.04", "", "Late", "r=0.35", "p=0.27"), pch=c(NA, NA, NA, NA, NA, NA, NA), lty=c(1, NA, NA, NA, 1, NA, NA), lwd=c(2,NA, NA, NA, 2, NA, NA), col=c("Dark Blue", NA, NA, NA, "Light Blue", NA, NA), cex=2)

fit2_late<-lm(NegChange$PercAngryEnd25_NonFace~NegChange$LogPostCort)

plot(NegChange[which(NegChange$Dx =="22q"),"LogPreCort"], NegChange[which(NegChange$Dx =="22q"),"PercAngryBegin25_NonFace"],col="Dark Blue", pch=20, cex=2, ylim=c(0.1,.8), xlim=c(-.5,.4), main = "Early and Late Gaze in Reducers \nNonFace Areas", xlab="Log Cortisol", ylab="Ratio Angry Face Viewing")
points(NegChange[which(NegChange$Dx =="TD"),"LogPreCort"], NegChange[which(NegChange$Dx =="TD"),"PercAngryBegin25_NonFace"], col="DarkBlue", pch=17, cex=2)
abline(fit2, lwd=2, col="Dark Blue")

points(NegChange[which(NegChange$Dx =="22q"),"LogPreCort"], NegChange[which(NegChange$Dx =="22q"),"PercAngryEnd25_NonFace"], col="Light Blue", pch=17, cex=2)
points(NegChange[which(NegChange$Dx =="TD"),"LogPreCort"], NegChange[which(NegChange$Dx =="TD"),"PercAngryEnd25_NonFace"], col="Light Blue", pch=17, cex=2)
abline(fit2_late, lwd=2, col="Light Blue")


plot(c(1,10),c(1,10))
legend(2, 9, c("Early", "r=-0.78", "p=0.003", "", "Late", "r=-0.60", "p=0.04"), pch=c(NA, NA, NA, NA, NA, NA, NA), lty=c(1, NA, NA, NA, 1, NA, NA), lwd=c(2,NA, NA, NA, 2, NA, NA), col=c("Dark Blue", NA, NA, NA, "Light Blue", NA, NA), cex=2)