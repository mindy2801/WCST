
rm(list=ls())

# read the data file
dat = read.table("reorganized data.txt", header=T)

allSubjs = unique(dat$subjID)  # all subject IDs
N = length(allSubjs)      # number of subjects
T = table(dat$subjID)     # number of trials per subject
numIter = 100             # number of iterations to find global minimum values
numPars = 3               # number of parameters

#Creating 3*4*128 answer sheet from "correctdeck matrix"
correctdeckmat=read.table("correctdeck matrix.txt",header=1) 
matchstack=array(rep(0,(3*4*128)),dim=c(3,4,128)) 
for (i in 1:4) matchstack[,i,] <- array(as.numeric(correctdeckmat[,]==i),dim=c(3,128))

#Calculate baseline model predction: deck-probabilities(4) based on attention(3)
cattpredp <- matrix(NA, ncol=128, nrow=4) #calculated predicted probabilities of choosing each deck under baseline model
for(trial in 1:128){
  cattpredp[,trial]=t(c(1/3,1/3,1/3)%*%matchstack[,,trial])
}




#Powerrize function. Necessary for parameter 'd'
powerrize=function(pow,tempvec) (tempvec^pow)/sum(tempvec^pow)

#Shaping parameter bounds
stretchpars=function(opars) -log((param_up-param_low)/(opars-param_low)-1)	#opars=original pars
contractpars=function(spars) (param_up-param_low)/(exp(-spars)+1)+param_low		#spars=stretched pars


#Return 2LL under model of interest(rpd1)
rpd1_mle = function(param, tmp_data){

  r = param[1]
  p = param[2]
  d = param[3]
  T = nrow(tmp_data)  # number of trials


  #Calculate probability for each deck.
  subjattmat=matrix(nrow=3,ncol=T) #initialize subject attention to 3 dimensions
  curatt=rep((1/3),3) 
  subjattmat[,1] <- curatt
  
  predpmat=matrix(-1,nrow=4, ncol=T) #initialize probability of choosing each deck(4 choices)
  predpmat[,1] <- t(curatt%*%matchstack[,,1])

    for (t in 1:(T-1)){
    attmatchchoice=correctdeckmat[,t]== tmp_data$deck[t] #for each choice, indicate which dimension matches the choice
    curatt <- curatt*.9999997+.0000001 #rescale to prevent rounding errors
    
      if(tmp_data$outcome[t]==1){
        attsignal <- powerrize(1,(attmatchchoice*curatt))
        curatt <- (1-r)*curatt+r*attsignal
      } else{
        attsignal <- powerrize(1,((1-attmatchchoice)*curatt))
        curatt <- (1-p)*curatt+p*attsignal
      }
    
      subjattmat[,t+1] <- curatt
      predpmat[,t+1] <- t(powerrize(d,curatt)%*%matchstack[,,t+1])
    }
  
  choiceprob <- c() #initialize probability of selected choices
  for(i in 1:T){
  choiceprob <- c(choiceprob, predpmat[tmp_data$deck[i],i])
  }
  choiceprob <- choiceprob*.9998+.0001
  -2*sum(log(choiceprob[-1])) #remove the first trial

}

#Return 2LL under baseline model
base_mle=function(tmp_data){
  T <- nrow(tmp_data)
  choiceprob <- c()
  
  for(i in 1:T){
    choiceprob <- c(choiceprob, cattpredp[tmp_data$deck[i],i])
  }
  choiceprob <- choiceprob*.9998+.0001
  
  return(-2*sum(log(choiceprob[-1])))
  
}

# parameter bounds
param_low <- c(0, 0, 0); param_up <- c(1, 1, 5);  # lower and upper bounds of r, p, d
param_initial = cbind(runif(numIter)*param_up[1], runif(numIter)*param_up[2], runif(numIter)*param_up[3])

global_pars = data.frame(r=NULL, p=NULL, d=NULL, AIC=NULL, BIC=NULL, BaseAIC=NULL, BaseBIC=NULL, subjID=NULL, group=NULL)

# for each subject.. 
for (i in 1:N) {
  tmp_ID = allSubjs[i]  # tmpID. e.g., tmpID = 2
  tmp_data = subset(dat, subjID == tmp_ID)  # select only tmp_ID's data
  tmp_group = as.character(unique(tmp_data$group))
  global_estimate = 1e6   # a very large number so that it can be rejected with an initial MLE estimate
  global_mle = NULL
  
  for (iter in 1:numIter) {
    # fine MLE estimates of the current subject
    mle_rpd1 <- optim(param_initial[iter, ], rpd1_mle, method="Nelder-Mead", tmp_data = tmp_data)
    #print(mle_ra_prospect$value) # for debugging, show outputs
    
    # Replace the results if the latest optimization yields better result
    if(mle_rpd1$value < global_estimate) {
      global_estimate <- mle_rpd1$value
      global_mle = mle_rpd1
    }
  }
  
  
  
  global_pars[i, "r"] = global_mle$par[1]
  global_pars[i, "p"] = global_mle$par[2]
  global_pars[i, "d"] = global_mle$par[3]
  global_pars[i, "AIC"] = global_mle$value + 2*numPars 
  global_pars[i, "BIC"] = global_mle$value + numPars*log(T[i]-1)
  global_pars[i, "Base_AIC"] = base_mle(tmp_data) + 2*numPars #HAVE TO CHECK FORMULA. What is the number of parameters of the baseline model?
  global_pars[i, "Base_BIC"] = base_mle(tmp_data) + numPars*log(T[i]-1)
  global_pars[i, "subjID"] = tmp_ID
  global_pars[i, "group"] = tmp_group
  
  
  cat("End of modeling subject ID =", tmp_ID, "\n")
}





print(global_pars)

sum_AIC = sum(global_pars$AIC)
print(sum_AIC)

sum_BIC = sum(global_pars$BIC)
print(sum_BIC)







###NEED TO CHECK
install.packages("numDeriv")
library(numDeriv)
?grad
grad(rpd1_mle, x=param_initial[10,], tmp_data=tmp_data)
grad2 <- function(param, tmp_data, fn=rpd1_mle) grad(x=param, func=fn, tmp_data=tmp_data)
optim(param_initial[1, ], fn=rpd1_mle, gr=grad2, method="L-BFGS-B",lower=param_low, upper=param_up, tmp_data=tmp_data) 
#I don't think I can use BGFS method here because likelihood is not derived from a continuous function.

#Error in convergence with L-BFGS-B
#Some problems in gradient of the MLE function?


#  $convergence
#[1] 0

#$message
#[1] "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"

#Error:
#  L-BFGS-B needs finite values of 'fn'
#non-finite finite-difference value 


###




