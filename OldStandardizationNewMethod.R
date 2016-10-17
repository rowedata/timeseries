library(Mcomp)
library(pROC)

rm(list = ls())# remove all objects

################  PART 1  ######################################
#put data in matrix form, where rows = observations, cols = variables 

X <- read.csv("12micro.csv", sep = ",", header = F)
p <- dim(X)
cols <- p[2]
row <- p[1]

########################################
#taking care of unequal length of the time series

rows <- matrix(, nrow = 1, ncol = cols) #calculating length of each of 474 series 

for (k in 1:cols){                      #and storing it as array, 
  rows [k] = row - sum(is.na(X[,k]) == 1)#can be referenced using column number
}

#####################################################################################
#simple exponential smoothing
P = 18               # 18 1-step ahead forecasts          
alpha = 0.1          # smoothing constant
fhats <- matrix(, nrow = P, ncol = cols) #matrix of forecasts

for (k in 1:cols){
  N = rows[k]- P                # the number of observations for twhich there is no AES prediction
                    
  prev = mean(X[1 : N-1, k]) 
  r = 1                         #create new index for predictions (1-18)
  
  for (i in (N + 1):rows[k]){
    error = X[i - 1, k] - prev  
    f  = prev + alpha * error # forecasted value for i+1-th observation,
    fhats[r,k] = f
    prev = f                  #updates for the next iteration
    r = r + 1
  }
  
}
 
####################################################################################
#adaptive exponential smoothing 


#preallocate matrix of predictions
fhat <- matrix(, nrow = P  , ncol = cols) #matrix of predictions
ds   <- matrix(, nrow = row, ncol = cols) # matrix of differences

for (k in 1 : cols){ 
  r = 1                         #counter for new predictions
  N = rows[k]- P                # the number of observations for twhich there is no AES prediction
  R = N + 1                     #number of first predicted observation 
      
      for (i in 2 : rows[k]) {            #create the differences
       ds[i, k] <-  X[i, k] - X[(i-1), k] 
      }
  
  stddev <- sd(ds[2 : N - 1], k) # std deviation of the differences
  
  prev = mean(X[1 : N - 1, k])   # starting value for mean 
  
  for (i in (N + 1) : rows[k]){     # i is the observation for which prediction is made
    
    error = X[i - 1, k] - prev # Xn -Xn-1 bar
    
    ratio = abs(error / stddev)  #ratio gives us idea if difference is too high/low
        
    if (ratio < 0.1) {     # "taylor's" part, weights and cutoffs picked by hand
          alpha = 0.05
        } else if (ratio < 1) {
          alpha = 0.1
        } else
          alpha = 0.2
    
    fhat[r,k] = prev + alpha * error # forecast with adjusted alpha
    
    prev = fhat [r, k]            #"base" for the next forecase - mean of all observations up to current period
    stddev = sd(ds[2 : i], k)           # std deviation is updated, all available observations are utilized
    r = r + 1
    
  }   
  
} 


##plot series against forecast

k = 70                  # pick a time series 
N = rows[k] - P + 1    # number of first predicted observation

v = X[(N):(rows[k]),k] # subsetting actuals and plotting them
plot(v, type = "l",xlim = c(1,P), xlab = k, ylab ="Value", sub ="simple - blue circles, new adaptive - red triangles" )

lines(fhat[,k], col = "red", type = "o", pch = 24) #add adaptive
lines(fhats[,k], col = "blue", type = "o")         #add simple

g = 5
cbind(X[(rows[g] - P + 1):rows[g], g], fhats[,g])

################ Accuracy measures for 2 methods#################################################
#######  ME   RMSE   MAE   MPE  MAPE  ##################

j1 <- matrix (, nrow = cols, ncol = 5) #simple
j2 <- matrix (, nrow = cols, ncol = 5) #adaptive

for (k in 1:cols){
  rt1 <- accuracy(fhats[,k], X[(rows[k] - P + 1):rows[k], k])
  rt2 <- accuracy(fhat[,k] , X[(rows[k] - P + 1):rows[k], k])
 
  j1[k,1:5] = rt1
  j2[k,1:5] = rt2
 
}

mean(j1[1:474, 2]) #simple RMSE
mean(j2[1:474, 2]) #better adaptive RMSE

mean(j1[1:474, 5]) #simple MAPE
mean(j2[1:474, 5]) #same adaptive MAPE


######Symmetric mean absolute precentage error (SMAPE) calculation##################
sums = 0
suma = 0


for (k in 1:cols){ 
    N = rows[k] - P
  
    for (i in 1:P){
      
      num = 2 * abs(fhats[i, k] - X[(N + i), k])
      den = fhats[i, k] + X[(N + i), k]
      j = num  / den
      sums = sums + j
      
      anum = 2 * abs(fhat[i, k] - X[(N + i), k])
      aden = fhat[i, k] + X[(N + i), k]
      p = anum  / aden 
      suma = suma + p
    } 
  
}

round(100 * sums /(cols*P),3) #symmetric MAPE for simple
round(100 * suma /(cols*P),3) # same symmetric MAPE adaptive



################  PART 2  ######################################
######### GENERATION OF RESPONSES ##############only last 18 obs are predicted

delt <- matrix(, nrow = row, ncol = cols) #create all possible deltas 

for (k in 1:cols){
  for (i in 2:rows[k]){
    delt[i, k] <- X[i, k]- X[i - 1, k] #current -previous, rows 2 till the end
  } 
}

#standardize actuals (delt)

## preallocate a matrix for standardized deltas 
zd  <- matrix(, nrow = P, ncol = cols)  # matrix of deltas
m   <- matrix(, nrow = 1, ncol = cols)  # means for each column
std <- matrix(, nrow = 1, ncol = cols)  # stds for each column
zrt <- matrix(nrow = 1, ncol= cols) #standardized delts for one period before forecast(used for RTM)


for (k in 1:cols){  
  
  end <- rows[k] - P #last observation used in prediction
  m[1,k] <- mean(delt[2 : end,k]) #one mean for all deltas in one column
  std[1,k] <- sd(delt[2 : end,k]) #one std for all deltas each series

      for (i in 1 : P){
        zd[i,k] <- (delt[(end + i), k] - m[1, k]) / (std[1,k])
      } 
  
  zrt[1,k] <- (delt[end,k]- m[1,k]) /std[1,k] #used for rtm
}




## create the gold standard

quantile = c(.95)
threshold = -qnorm(quantile)

#preallocate matrix of responses
response <-  matrix (0, nrow = P, ncol = cols ) # matrix of responses of size 40 by 474

# get responses
for (k in 1:cols){
  
  for (i in 1:P){
    if (zd[i,k] < threshold) {
      response[i,k] = 1          #signal= 1 if below negative threshold, 
    } 
  }
  
}

sum(response[1,], na.rm = T) #number of positives by row without RTM

#check regression to the mean (RTM)

for (k in 1:cols){
    if ( zrt[1,k] > abs(threshold) ) {
    
    response[1, k] = 0
  }
}

sum(response[1,], na.rm = T) #number of positives by row with RTM


########################################################
#OLD STANDARDIZATION PROCEDURE


adelt <- matrix(, nrow = P, ncol = cols) #adaptive forecast delta
sdelt <- matrix(, nrow = P, ncol = cols) #simple forecast delta

for (k in 1:cols){
  
  N  = rows[k]- P-1  #row number one before the first forecast
  
  for (i in 1:P){
    
    adelt[i, k] <- (fhat[i, k]  - X[N + i, k])    #forecasted change AES 
    sdelt[i, k] <- (fhats[i, k]  - X[N + i, k])   #forecasted change SES
  }
}


## standardized deltas for forecasts
za <- matrix(, nrow = P, ncol = cols) #adaptive
zs <- matrix(, nrow = P, ncol = cols) #simple

for (k in 1:cols){ 
  end <- rows[k] - P
  
  for (i in 1:P){
    
    m  <- mean(delt[2 : end,k]) #here we're using proportionality 
    std<- sd(delt[2 : end,k]) #i.e. deltas are used to standardize a/sdeltas
   
    za[i, k] <- (adelt[i, k] - m)/ std
    zs[i, k] <- (sdelt[i, k] - m)/ std
    
  }
}


######################################################################################
################  PART 3  ######################################
############### ROC CALCULATION ###########################################################

## let's check the column means and std. devs. 
round(colMeans(zd, na.rm = T), 3)      # seems near 0
round(apply(zd, 2, sd, na.rm = T) , 3)  # should be all 1's or close

k = 67 ## okay so now we want to compare the actual 1s and 0s to the forecasts ## response cols should now have a 0 or 1 , or NA if no data there

j <- roc(na.omit(response[,k]), na.omit(za[,k]),plot = T)# this would be an ROC curve for 1 time series
roc.test(na.omit(response[,k]), na.omit(za[,k]),na.omit(zs[,k])) ##  test whether za is better than zs

cbind(response[,7], round(zd[,7],2)) # see corresspondence between responses and scores

################################################################################
#combine it all into one vector and calculate ROC

t = 474 #calculate the total number of responses from dataset

#combine all responses and all scores for SES and AES in 3 separate vectors: RESPone,SESone, AESone

RESP1 <- response[1,]
SES1 <- zs[1,]
AES1 <- za[1,]
TES1 <- zt[1,]



AES <- roc(RESP1, AES1, plot = T, main = "Adaptive ES", panel.first = grid(lwd = 1))
SES <- roc(RESP1, SES1, plot = T, main = "Simple ES", panel.first = grid(lwd = 1))

#plot two ROC curves together
plot(AES,  col = "red", legacy.axes = TRUE, smooth=FALSE)
plot(SES, add = TRUE, smooth=FALSE)
grid (5,10, lty = 8) 
legend("bottomright", col=c("black","red" ),legend =c("simple ES","adaptive ES"), ncol = 1, lty =1 )

#calculate PAUCs for all made predictions
pauc1 <- auc(RESP1, SES1, partial.auc = c(1, 0.8)) #simple
pauc2 <- auc(RESP1, AES1, partial.auc = c(1, 0.8)) #adaptive


round(pauc1,3) #simple pauc
round(pauc2,3) #adaptive pauc

SES       
AES

#significance testing 
temp=roc.test(RESP1,SES1,AES1,method="bootstrap", boot.n=1000,partial.auc=c(1, 0.8))
temp

#partial ROC curve for SES

