library(Mcomp)
library(pROC)

quantilevect=c(.99,.95,.90)
horizonvect=c(1,2,4)
for (v in 1:3){
for (w in 1:3){
quantile=quantilevect[v]
threshold=-qnorm(quantile)  ## gold standard cutoff
horizon=horizonvect[w]  ## forecast horizon for analysis
series=474  ## number of time series
numf=18 ## number of forecasts (in M3 = 18)
numa=108 ## max number of past actuals (in M3 = 108)

mm<-subset(M3,12,"micro")  ##1402 to 1875

##create a matrix for the forecasts
holdout=rep(-1,numf*series)
attr(holdout,"dim")<-c(numf,series)
##index=1402
##paste("mm$N",index,"$xx[",j,"]",sep="")
holdoutdiff=holdout

## Put all forecast horizons into a matrix
for (i in 1:series){
temp=mm[i]
temp=temp$N
xx=temp$xx
for (j in 1:numf){
holdout[j,i]=xx[j]
}}



### create a matrix for the actuals

actuals=rep(NA,numa*series)
attr(actuals,"dim")<-c(numa,series)
actualsdiff=actuals

### Put all actuals into a matrix and line up vectors on the bottom. empties on top.
for (i in 1:series){
temp=mm[i]
temp=temp$N
x=temp$x
for (j in 1:length(x)){
actuals[numa+1-j,i]=x[length(x)+1-j]
}}

##fill in differences (Actual horizon i minus last actual)
for (i in 1:numf){
holdoutdiff[i,]=holdout[i,]-actuals[numa,]
}


## fill in the matrix of past deltas
actualsdiff=diff(actuals)

## create a matrix for standardized differences
holdoutdiffstd=holdoutdiff
stds=rep(0,series)
for (i in 1:series){
stds[i]=sd(actualsdiff[,i],na.rm=T)
}

for (i in 1:numf){
holdoutdiffstd[i,]=(holdoutdiffstd[i,]-as.vector(colMeans(actualsdiff,na.rm=T)))/as.vector(stds)
}

## create the gold standard

gs=rep(0,series)
for (i in 1:series){
if (holdoutdiffstd[horizon,i]<threshold) gs[i]=1
}
gsraw<-sum(gs)


##regression to the mean
##cancel a gold standard if previous direction was opposite


rgttm=(actualsdiff[numa-1,]-as.vector(colMeans(actualsdiff,na.rm=T)))/as.vector(stds)

for (i in 1:series){
if (rgttm[i]>abs(threshold)) gs[i]=0
}

gsrtm<-sum(gs)


###get the forecastfiles
forecast=read.csv("bigforecast2.csv",sep=",",header=F)

######################################
#################################
#############CHECK#######################
###############################
##forecast2=read.table("Rforecasth1.csv",sep=",")
##forecasthorizon2=matrix(-1,nrow=24,ncol=474)

##for (k in 1:474){
##for (i in 1:24){
##forecasthorizon2[i,k]=forecast2[i,k]
##}}

##forecasthorizon2=forecasthorizon2[-23,]

##test=forecasthorizon-forecasthorizon2
#################################
###################################

##make a matrix of only that forecast horizon
j=horizon
forecasthorizon=matrix(-1,nrow=24,ncol=474)

for (k in 1:474){
for (i in 1:24){
forecasthorizon[i,k]=forecast[j+(i-1)*18,k]
}}


temp=as.vector(actuals[numa,])
temp2=as.vector(colMeans(actualsdiff,na.rm=T))
temp3=as.vector(stds)

##standardize the forecast matrix of horizon j
for (i in 1:24){
forecasthorizon[i,]=forecasthorizon[i,]-temp
forecasthorizon[i,]=(forecasthorizon[i,]-temp2)/temp3
}

## j is horizon  and i is forecasting method
##compute pauc of .2

##create combination method forecasts using  Forecast Pro, Automattann, and Theta
## these are rows 21, 18, and 22
### put in place of AAM1, AAM2, and rule based forecasting
## these are rows 15, 16, and 1




#########################
#########################
####RESULT PORTION#######
#########################
#########################



names=c("RULE BASED F","AUTOBOX2","AUTOBOX3","AUTOBOX1","ARARMA","BJ AUTOMATIC","COM S-H-D","DAMPEN","FLORES PEARCE 1","FLORES PEARCE 2","PP AUTOCAST","HOLT","FORECAST X","NAÏVE 2","AAM1","AAM2","ROBUST TREND","AUTOMATTANN","SINGLE","SMART FCS","FORECAST PRO","THETA","THETA SM","WINTER")
paucs=rep(0,24)

for (i in 1:24){
paucs[i]=auc(gs,forecasthorizon[i,],partial.auc=c(1,.8))

print(c(names[i],paucs[i]))
}


paucs=round(paucs,digits=4)

for (i in 1:24){
print(c(names[i],paucs[i]))
}

ranks=matrix(nrow=24,ncol=2)
ranks[,1]=names
ranks[,2]=paucs
ranks[,2]=as.numeric(ranks[,2])
ranks=ranks[order(ranks[,2],decreasing=TRUE),] 
ind=sort(paucs,decreasing=T,index.return = TRUE)
ix=ind$ix

temp2=rep(-1,24)
temp2[1]=NA
for (j in 2:24){
temp=roc.test(gs,forecasthorizon[ix[1],],forecasthorizon[ix[j],],method="bootstrap", boot.n=1000,partial.auc=c(1, 0.8))
temp2[j]=temp$p.value
}
temp3=round(temp2,digits=4)
ranks=cbind(ranks,temp3,temp3/2)
ranks[1,3]="p-value"
ranks[1,4]="one-sided"

write.csv(ranks, paste(file="pauc",quantile,"horizon",horizon,".csv",sep=""), row.names = FALSE)




#########################
###Kendall's Tau#########
######################
library(Kendall)
complex=c(2,3,4.3,5.3,6.2,7.5,8,8,11,11.3,12,14.3,14.3,15.2,16.8,16.8,16.8,17.5,18.8,21.5)
namesc=c("SINGLE","HOLT","ROBUST TREND","WINTER","DAMPEN",
"PP AUTOCAST","THETA SM","COM S-H-D","THETA","BJ AUTOMATIC",
"AUTOBOX1","AUTOBOX2","AUTOBOX3","ARARMA","SMART FCS",
"FLORES PEARCE 2","FLORES PEARCE 1","FORECAST PRO",
"FORECAST X","AUTOMATTANN")
kcomplex=matrix(nrow=3,ncol=20)
kcomplex[1,]=namesc
kcomplex[2,]=complex
for (i in 1:20){
kcomplex[3,i]=ranks[which(ranks[,1]==namesc[i]),2]
}
print(summary(Kendall(kcomplex[2,],kcomplex[3,])))
#####################
#####################

print(gsraw)
print(gsrtm)




}}

