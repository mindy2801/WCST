
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

#Powerrize function
powerrize=function(pow,tempvec) (tempvec^pow)/sum(tempvec^pow)

#Return 2LL
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
  
  choiceprob <- c() #initialize probability of actual choices
  for(i in 1:T){
  choiceprob <- c(choiceprob, predpmat[tmp_data$deck[i],i])
  }
  -2*sum(log(choiceprob[-1])) #remove the first trial
  
  

}

# parameter bounds
param_low <- c(0, 0, 0); param_up <- c(1, 1, 5);  # lower and upper bounds of r, p, d

# randomly generate initial values for numIter iterations. 
# param_initial --> a matrix with numIter*3 dimension. Use this for all subjects.
param_initial = cbind(runif(numIter)*param_up[1], runif(numIter)*param_up[2], runif(numIter)*param_up[3])

global_pars = data.frame(r=NULL, p=NULL, d=NULL, AIC=NULL, BIC=NULL, subjID=NULL, group=NULL)

# for each subject.. 
for (i in 1:N) {
  tmp_ID = allSubjs[i]  # tmpID. e.g., tmpID = 2
  tmp_data = subset(dat, subjID == tmp_ID)  # select only tmp_ID's data
  tmp_group = as.character(unique(tmp_data$group))
  global_estimate = 1e6   # a very large number so that it can be rejected with an initial MLE estimate
  global_mle = NULL
  
  for (iter in 1:numIter) {
    # fine MLE estimates of the current subject
    mle_rpd1 <- optim(param_initial[iter, ], rpd1_mle, method="L-BFGS-B", 
                             lower=param_low, upper=param_up, tmp_data = tmp_data)
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
  global_pars[i, "AIC"] = global_mle$value + 2*numPars #HAVE TO CHECK FORMULA
  global_pars[i, "BIC"] = global_mle$value + numPars*log(T-1)
  global_pars[i, "subjID"] = tmp_ID
  global_pars[i, "group"] = tmp_group
  
  cat("End of modeling subject ID =", tmp_ID, "\n")
}

print(global_pars)

sum_AIC = sum(global_pars$AIC)
print(sum_AIC)

sum_BIC = sum(global_pars$BIC)
print(sum_BIC)







