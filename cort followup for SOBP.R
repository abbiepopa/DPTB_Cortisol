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
 
 ###Look at low General ABAS vs. high General ABAS
 AllDataCortAnalysis$LogicABASGACImpaired<-(AllDataCortAnalysis$ABAS_GAC<76)
 ###Look at low Total Spence vs high Total Spence
 AllDataCortAnalysis$Spnc_Eleveated<-(AllDataCortAnalysis$Spnc_Ttl_Parent>59)
 havecort<-subset(AllDataCortAnalysis, !is.na(AllDataCortAnalysis$PreCORT))

 table(havecort$LogicABASGACImpaired, havecort$Dx)

 table(havecort$LogicABASGACImpaired, havecort$Cluster)
       
 table(havecort$Spnc_Eleveated, havecort$Dx)
       
 table(havecort$Spnc_Eleveated, havecort$Cluster)
       
 table(havecort$Dx)
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

 havecort<-subset(AllDataCortAnalysis, !is.na(AllDataCortAnalysis$PreCORT))

#Angry ABAS

HighABAS_AngryNonFace<-cor.test(havecort[which(!havecort$LogicABASGACImpaired),]$CortDelta, havecort[which(!havecort$LogicABASGACImpaired),]$PercAngryNonFace)
LowABAS_AngryNonFace<-cor.test(havecort[which(havecort$LogicABASGACImpaired),]$CortDelta, havecort[which(havecort$LogicABASGACImpaired),]$PercAngryNonFace)

zdiff_ABAS_AngryNonFace<-(HighABAS_AngryNonFace$estimate[[1]] - LowABAS_AngryNonFace$estimate[[1]])/sqrt((1/30)+(1/27))


HighABAS_AngryAngry<-cor.test(havecort[which(!havecort$LogicABASGACImpaired),]$CortDelta, havecort[which(!havecort$LogicABASGACImpaired),]$PercAngry)
LowABAS_AngryAngry<-cor.test(havecort[which(havecort$LogicABASGACImpaired),]$CortDelta, havecort[which(havecort$LogicABASGACImpaired),]$PercAngry)

zdiff_ABAS_AngryAngry<-(HighABAS_AngryAngry$estimate[[1]] - LowABAS_AngryAngry$estimate[[1]])/sqrt((1/30)+(1/27))

#Happy ABAS

HighABAS_HappyNonFace<-cor.test(havecort[which(!havecort$LogicABASGACImpaired),]$CortDelta, havecort[which(!havecort$LogicABASGACImpaired),]$PercHappyNonFace)
LowABAS_HappyNonFace<-cor.test(havecort[which(havecort$LogicABASGACImpaired),]$CortDelta, havecort[which(havecort$LogicABASGACImpaired),]$PercHappyNonFace)

zdiff_ABAS_HappyNonFace<-(HighABAS_HappyNonFace$estimate[[1]] - LowABAS_HappyNonFace$estimate[[1]])/sqrt((1/30)+(1/27))


HighABAS_HappyHappy<-cor.test(havecort[which(!havecort$LogicABASGACImpaired),]$CortDelta, havecort[which(!havecort$LogicABASGACImpaired),]$PercHappy)
LowABAS_HappyHappy<-cor.test(havecort[which(havecort$LogicABASGACImpaired),]$CortDelta, havecort[which(havecort$LogicABASGACImpaired),]$PercHappy)

zdiff_ABAS_HappyHappy<-(HighABAS_HappyHappy$estimate[[1]] - LowABAS_HappyHappy$estimate[[1]])/sqrt((1/30)+(1/27))

#Angry Spence

LowSpence_AngryNonFace<-cor.test(havecort[which(!havecort$Spnc_Eleveated),]$CortDelta, havecort[which(!havecort$Spnc_Eleveated),]$PercAngryNonFace)
HighSpence_AngryNonFace<-cor.test(havecort[which(havecort$Spnc_Eleveated),]$CortDelta, havecort[which(havecort$Spnc_Eleveated),]$PercAngryNonFace)

zdiff_Spence_AngryNonFace<-(HighSpence_AngryNonFace$estimate[[1]] - LowSpence_AngryNonFace$estimate[[1]])/sqrt((1/29)+(1/31))


LowSpence_AngryNeutral<-cor.test(havecort[which(!havecort$Spnc_Eleveated),]$CortDelta, havecort[which(!havecort$Spnc_Eleveated),]$PercAngryNeutral)
HighSpence_AngryNeutral<-cor.test(havecort[which(havecort$Spnc_Eleveated),]$CortDelta, havecort[which(havecort$Spnc_Eleveated),]$PercAngryNeutral)

zdiff_Spence_AngryNeutral<-(HighABAS_AngryAngry$estimate[[1]] - LowABAS_AngryAngry$estimate[[1]])/sqrt((1/29)+(1/31))

#Happy Spence

LowSpence_HappyNonFace<-cor.test(havecort[which(!havecort$Spnc_Eleveated),]$CortDelta, havecort[which(!havecort$Spnc_Eleveated),]$PercHappyNonFace)
HighSpence_HappyNonFace<-cor.test(havecort[which(havecort$Spnc_Eleveated),]$CortDelta, havecort[which(havecort$Spnc_Eleveated),]$PercHappyNonFace)

zdiff_Spence_HappyNonFace<-(HighABAS_HappyNonFace$estimate[[1]] - LowABAS_HappyNonFace$estimate[[1]])/sqrt((1/29)+(1/31))


LowSpence_HappyHappy<-cor.test(havecort[which(!havecort$Spnc_Eleveated),]$CortDelta, havecort[which(!havecort$Spnc_Eleveated),]$PercHappy)
HighSpence_HappyHappy<-cor.test(havecort[which(havecort$Spnc_Eleveated),]$CortDelta, havecort[which(havecort$Spnc_Eleveated),]$PercHappy)

zdiff_Spence_HappyHappy<-(HighABAS_HappyHappy$estimate[[1]] - LowABAS_HappyHappy$estimate[[1]])/sqrt((1/29)+(1/31))


###PreCort
#Angry Changes, overall faces

NegChange_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryNonFace)

NoChange_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==0),]$LogPreCort, havecort[which(havecort$CortDeltaBin==0),]$PercAngryNonFace)

PosChange_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==1),]$PercAngryNonFace)

zdiff_NegNo_AngryNonFace<-(NegChange_AngryNonFace$estimate[[1]]-NoChange_AngryNonFace$estimate[[1]])/sqrt((1/13)+(1/32))

zdiff_NegPos_AngryNonFace<-(NegChange_AngryNonFace$estimate[[1]]-PosChange_AngryNonFace$estimate[[1]])/sqrt((1/13)+(1/12))

zdiff_NoPos_AngryNonFace<-(NoChange_AngryNonFace$estimate[[1]]-PosChange_AngryNonFace$estimate[[1]])/sqrt((1/32)+(1/12))

#Happy Changes, overall faces
#nonface

NegChange_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyNonFace)

NoChange_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==0),]$LogPreCort, havecort[which(havecort$CortDeltaBin==0),]$PercHappyNonFace)

PosChange_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==1),]$PercHappyNonFace)

zdiff_NegNo_HappyNonFace<-(NegChange_HappyNonFace$estimate[[1]]-NoChange_HappyNonFace$estimate[[1]])/sqrt((1/13)+(1/32))

zdiff_NegPos_HappyNonFace<-(NegChange_HappyNonFace$estimate[[1]]-PosChange_HappyNonFace$estimate[[1]])/sqrt((1/13)+(1/12))

zdiff_NoPos_HappyNonFace<-(NoChange_HappyNonFace$estimate[[1]]-PosChange_HappyNonFace$estimate[[1]])/sqrt((1/32)+(1/12))

#happy

NegChange_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappy)

NoChange_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==0),]$LogPreCort, havecort[which(havecort$CortDeltaBin==0),]$PercHappy)

PosChange_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==1),]$PercHappy)

zdiff_NegNo_HappyHappy<-(NegChange_HappyHappy$estimate[[1]]-NoChange_HappyHappy$estimate[[1]])/sqrt((1/13)+(1/32))

zdiff_NegPos_HappyHappy<-(NegChange_HappyHappy$estimate[[1]]-PosChange_HappyHappy$estimate[[1]])/sqrt((1/13)+(1/12))

zdiff_NoPos_HappyHappy<-(NoChange_HappyHappy$estimate[[1]]-PosChange_HappyHappy$estimate[[1]])/sqrt((1/32)+(1/12))

PreValues<-list(NegChange_AngryNonFace, NoChange_AngryNonFace, PosChange_AngryNonFace, zdiff_NegNo_AngryNonFace, zdiff_NegPos_AngryNonFace, zdiff_NoPos_AngryNonFace, NegChange_HappyNonFace, NoChange_HappyNonFace, PosChange_HappyNonFace, zdiff_NegNo_HappyNonFace, zdiff_NegPos_HappyNonFace, zdiff_NoPos_HappyNonFace, NegChange_HappyHappy, NoChange_HappyHappy, PosChange_HappyHappy, zdiff_NegNo_HappyHappy, zdiff_NegPos_HappyHappy, zdiff_NoPos_HappyHappy)

###PostCort
#Angry Changes, overall faces

NegChange_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryNonFace)

NoChange_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==0),]$LogPostCort, havecort[which(havecort$CortDeltaBin==0),]$PercAngryNonFace)

PosChange_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==1),]$PercAngryNonFace)

zdiff_NegNo_AngryNonFace<-(NegChange_AngryNonFace$estimate[[1]]-NoChange_AngryNonFace$estimate[[1]])/sqrt((1/13)+(1/32))

zdiff_NegPos_AngryNonFace<-(NegChange_AngryNonFace$estimate[[1]]-PosChange_AngryNonFace$estimate[[1]])/sqrt((1/13)+(1/12))

zdiff_NoPos_AngryNonFace<-(NoChange_AngryNonFace$estimate[[1]]-PosChange_AngryNonFace$estimate[[1]])/sqrt((1/32)+(1/12))

#Happy Changes, overall faces
#nonface

NegChange_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyNonFace)

NoChange_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==0),]$LogPostCort, havecort[which(havecort$CortDeltaBin==0),]$PercHappyNonFace)

PosChange_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==1),]$PercHappyNonFace)

zdiff_NegNo_HappyNonFace<-(NegChange_HappyNonFace$estimate[[1]]-NoChange_HappyNonFace$estimate[[1]])/sqrt((1/13)+(1/32))

zdiff_NegPos_HappyNonFace<-(NegChange_HappyNonFace$estimate[[1]]-PosChange_HappyNonFace$estimate[[1]])/sqrt((1/13)+(1/12))

zdiff_NoPos_HappyNonFace<-(NoChange_HappyNonFace$estimate[[1]]-PosChange_HappyNonFace$estimate[[1]])/sqrt((1/32)+(1/12))

#happy

NegChange_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappy)

NoChange_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==0),]$LogPostCort, havecort[which(havecort$CortDeltaBin==0),]$PercHappy)

PosChange_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==1),]$PercHappy)

zdiff_NegNo_HappyHappy<-(NegChange_HappyHappy$estimate[[1]]-NoChange_HappyHappy$estimate[[1]])/sqrt((1/13)+(1/32))

zdiff_NegPos_HappyHappy<-(NegChange_HappyHappy$estimate[[1]]-PosChange_HappyHappy$estimate[[1]])/sqrt((1/13)+(1/12))

zdiff_NoPos_HappyHappy<-(NoChange_HappyHappy$estimate[[1]]-PosChange_HappyHappy$estimate[[1]])/sqrt((1/32)+(1/12))

PostValues<-list(NegChange_AngryNonFace, NoChange_AngryNonFace, PosChange_AngryNonFace, zdiff_NegNo_AngryNonFace, zdiff_NegPos_AngryNonFace, zdiff_NoPos_AngryNonFace, NegChange_HappyNonFace, NoChange_HappyNonFace, PosChange_HappyNonFace, zdiff_NegNo_HappyNonFace, zdiff_NegPos_HappyNonFace, zdiff_NoPos_HappyNonFace, NegChange_HappyHappy, NoChange_HappyHappy, PosChange_HappyHappy, zdiff_NegNo_HappyHappy, zdiff_NegPos_HappyHappy, zdiff_NoPos_HappyHappy)

###PreCort
###Negative Change, but comparing beginning, overall, and end
#Angry
Begin_AngryNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryBegin25_Neutral)
Overall_AngryNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryNeutral)
End_AngryNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryEnd25_Neutral)

zdiff_BeginOverall_AngryNeutral<-(Begin_AngryNeutral$estimate[[1]]-Overall_AngryNeutral$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_AngryNeutral<-(Begin_AngryNeutral$estimate[[1]]-End_AngryNeutral$estimate[[1]])/sqrt((1/13)+(1/13))

Begin_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryBegin25_NonFace)
Overall_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryNonFace)
End_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryEnd25_NonFace)

zdiff_BeginOverall_AngryNonFace<-(Begin_AngryNonFace$estimate[[1]]-Overall_AngryNonFace$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_AngryNonFace<-(Begin_AngryNonFace$estimate[[1]]-End_AngryNonFace$estimate[[1]])/sqrt((1/13)+(1/13))

Begin_AngryAngry<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryBegin25_Angry)
Overall_AngryAngry<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngry)
End_AngryAngry<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryEnd25_Angry)

zdiff_BeginOverall_AngryAngry<-(Begin_AngryAngry$estimate[[1]]-Overall_AngryAngry$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_AngryAngry<-(Begin_AngryAngry$estimate[[1]]-End_AngryAngry$estimate[[1]])/sqrt((1/13)+(1/13))

#Happy
Begin_HappyNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyBegin25_Neutral)
Overall_HappyNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyNeutral)
End_HappyNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyEnd25_Neutral)

zdiff_BeginOverall_HappyNeutral<-(Begin_HappyNeutral$estimate[[1]]-Overall_HappyNeutral$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_HappyNeutral<-(Begin_HappyNeutral$estimate[[1]]-End_HappyNeutral$estimate[[1]])/sqrt((1/13)+(1/13))

Begin_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyBegin25_NonFace)
Overall_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyNonFace)
End_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyEnd25_NonFace)

zdiff_BeginOverall_HappyNonFace<-(Begin_HappyNonFace$estimate[[1]]-Overall_HappyNonFace$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_HappyNonFace<-(Begin_HappyNonFace$estimate[[1]]-End_HappyNonFace$estimate[[1]])/sqrt((1/13)+(1/13))

Begin_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyBegin25_Happy)
Overall_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappy)
End_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPreCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyEnd25_Happy)

zdiff_BeginOverall_HappyHappy<-(Begin_HappyHappy$estimate[[1]]-Overall_HappyHappy$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_HappyHappy<-(Begin_HappyHappy$estimate[[1]]-End_HappyHappy$estimate[[1]])/sqrt((1/13)+(1/13))

preneglist<-list(Begin_AngryNeutral,Overall_AngryNeutral, End_AngryNeutral, zdiff_BeginOverall_AngryNeutral, zdiff_BeginEnd_AngryNeutral, Begin_AngryNonFace, Overall_AngryNonFace, End_AngryNonFace, zdiff_BeginOverall_AngryNonFace, zdiff_BeginEnd_AngryNonFace, Begin_AngryAngry, Overall_AngryAngry, End_AngryAngry, zdiff_BeginOverall_AngryAngry, zdiff_BeginEnd_AngryAngry, Begin_HappyNeutral, Overall_HappyNeutral, End_HappyNeutral, zdiff_BeginOverall_HappyNeutral, zdiff_BeginEnd_HappyNeutral, Begin_HappyNonFace, Overall_HappyNonFace, End_HappyNonFace, zdiff_BeginOverall_HappyNonFace, zdiff_BeginEnd_HappyNonFace, Begin_HappyHappy, Overall_HappyHappy, End_HappyHappy, zdiff_BeginOverall_HappyHappy, zdiff_BeginEnd_HappyHappy)

###PostCort
###Negative Change, but comparing beginning, overall, and end
#Angry
Begin_AngryNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryBegin25_Neutral)
Overall_AngryNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryNeutral)
End_AngryNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryEnd25_Neutral)

zdiff_BeginOverall_AngryNeutral<-(Begin_AngryNeutral$estimate[[1]]-Overall_AngryNeutral$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_AngryNeutral<-(Begin_AngryNeutral$estimate[[1]]-End_AngryNeutral$estimate[[1]])/sqrt((1/13)+(1/13))

Begin_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryBegin25_NonFace)
Overall_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryNonFace)
End_AngryNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryEnd25_NonFace)

zdiff_BeginOverall_AngryNonFace<-(Begin_AngryNonFace$estimate[[1]]-Overall_AngryNonFace$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_AngryNonFace<-(Begin_AngryNonFace$estimate[[1]]-End_AngryNonFace$estimate[[1]])/sqrt((1/13)+(1/13))

Begin_AngryAngry<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryBegin25_Angry)
Overall_AngryAngry<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngry)
End_AngryAngry<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryEnd25_Angry)

zdiff_BeginOverall_AngryAngry<-(Begin_AngryAngry$estimate[[1]]-Overall_AngryAngry$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_AngryAngry<-(Begin_AngryAngry$estimate[[1]]-End_AngryAngry$estimate[[1]])/sqrt((1/13)+(1/13))

#Happy
Begin_HappyNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyBegin25_Neutral)
Overall_HappyNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyNeutral)
End_HappyNeutral<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyEnd25_Neutral)

zdiff_BeginOverall_HappyNeutral<-(Begin_HappyNeutral$estimate[[1]]-Overall_HappyNeutral$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_HappyNeutral<-(Begin_HappyNeutral$estimate[[1]]-End_HappyNeutral$estimate[[1]])/sqrt((1/13)+(1/13))

Begin_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyBegin25_NonFace)
Overall_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyNonFace)
End_HappyNonFace<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyEnd25_NonFace)

zdiff_BeginOverall_HappyNonFace<-(Begin_HappyNonFace$estimate[[1]]-Overall_HappyNonFace$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_HappyNonFace<-(Begin_HappyNonFace$estimate[[1]]-End_HappyNonFace$estimate[[1]])/sqrt((1/13)+(1/13))

Begin_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyBegin25_Happy)
Overall_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappy)
End_HappyHappy<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercHappyEnd25_Happy)

zdiff_BeginOverall_HappyHappy<-(Begin_HappyHappy$estimate[[1]]-Overall_HappyHappy$estimate[[1]])/sqrt((1/13)+(1/13))
zdiff_BeginEnd_HappyHappy<-(Begin_HappyHappy$estimate[[1]]-End_HappyHappy$estimate[[1]])/sqrt((1/13)+(1/13))

postneglist<-list(Begin_AngryNeutral,Overall_AngryNeutral, End_AngryNeutral, zdiff_BeginOverall_AngryNeutral, zdiff_BeginEnd_AngryNeutral, Begin_AngryNonFace, Overall_AngryNonFace, End_AngryNonFace, zdiff_BeginOverall_AngryNonFace, zdiff_BeginEnd_AngryNonFace, Begin_AngryAngry, Overall_AngryAngry, End_AngryAngry, zdiff_BeginOverall_AngryAngry, zdiff_BeginEnd_AngryAngry, Begin_HappyNeutral, Overall_HappyNeutral, End_HappyNeutral, zdiff_BeginOverall_HappyNeutral, zdiff_BeginEnd_HappyNeutral, Begin_HappyNonFace, Overall_HappyNonFace, End_HappyNonFace, zdiff_BeginOverall_HappyNonFace, zdiff_BeginEnd_HappyNonFace, Begin_HappyHappy, Overall_HappyHappy, End_HappyHappy, zdiff_BeginOverall_HappyHappy, zdiff_BeginEnd_HappyHappy)

Neg_End_AngryAngry<-cor.test(havecort[which(havecort$CortDeltaBin==-1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==-1),]$PercAngryEnd25_Angry)
Pos_End_AngryAngry<-cor.test(havecort[which(havecort$CortDeltaBin==1),]$LogPostCort, havecort[which(havecort$CortDeltaBin==1),]$PercAngryEnd25_Angry)
zdiff_NegPos_AAEnd<-(Neg_End_AngryAngry$estimate[[1]]-Pos_End_AngryAngry[[1]])/sqrt((1/13)+(1/12))