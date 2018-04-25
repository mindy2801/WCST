
rm(list=ls())

# read the data file
dat = read.table("reorganized data.txt", header=T)

allSubjs = unique(dat$subjID)  # all subject IDs
N = length(allSubjs)      # number of subjects
T = table(dat$subjID)     # number of trials per subject
numIter = 100             # number of iterations to find global minimum values
numPars = 3               # number of parameters

WCST_mle = function(param, tmp_data)  {
  # param --> a vector with 3 initial parameter values. e.g., param = c(0.5, 1, 1)
  # tmp_data --> a single subject's data (140 trials), containing at least the following columns
  #              gain, loss, cert, gamble
  # 
  # evSafe: expected value of a certain (safe) option
  # evGamble: expected value of a risky option (gamble)
  # pGamble   # probability of choosing a gamble on each trial
  # free parameters: rho, tau, lambda
  
  r = param[1]
  p = param[2]
  d = param[3]
  T = nrow(tmp_data)  # number of trials
  sum_minusLL = 0  # sum of minus log likelihood. Initialize
  ###임시!
  #Calculate probability for each deck.
  vattpredpfun9=function(param,deck,outcome){
    
    curatt=rep((1/3),3)
    subjattmat=matrix(nrow=3,ncol=T) #subject attention to 3 dimensions
    predpmat=matrix(-1,nrow=4, ncol=T) #probability of choosing each deck(4 choices)
    subjattmat[,1]=curatt
    for (temptrial in 1:(T-1))
    {
      attmatchchoice=correctdeckmat[,temptrial]== tempchoices[,temptrial]
      doubleatt=curatt*.9999997+.0000001 #rescale to prevent rounding errors
      if(tempreinf[,temptrial]==1) 
      {
        attsignal=powerrize(i,(attmatchchoice*doubleatt))
        curatt=(1-r)*curatt+r*attsignal
      } else
      {
        attsignal=powerrize(i,((1-attmatchchoice)*doubleatt))
        curatt=(1-p)*curatt+p*attsignal
      }
      predpmat[,temptrial+1]=t(powerrize(d,curatt)%*%matchstack[,,temptrial+1])
    }
    temppars <- c(r,p,d)
    return(predpmat)
  }
  
  
  vattG2overarchfun=function(temppars,tempparbounds,tempchoices,tempreinf,predpfun) 
  {
    origpars=contractpars(temppars)
    templength=length(tempchoices)
    tempchoiceprob=matrix(c(tempchoices==1,tempchoices==2,tempchoices==3,tempchoices==4),nrow=4,ncol=templength,byrow=TRUE)*( predpfun(origpars,tempchoices,tempreinf))
    tempchoiceprob=colSums(tempchoiceprob)
    tempchoiceprob=tempchoiceprob[2:templength]*.9998+.0001		#removes trial 1 and rescales
    return(-2*sum(log(tempchoiceprob)))
  }
  
  ###임시!
  cert = tmp_data$cert
  gain = tmp_data$gain
  loss = tmp_data$loss
  gamble = tmp_data$gamble
  
  for (t in 1:T) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)
}

# parameter bounds
param_low <- c(0, 0, 0); param_up <- c(2, 10, 5);  # lower and upper bounds of ra_prospect model

# randomly generate initial values for numIter iterations. 
# param_initial --> a matrix with numIter*3 dimension. Use this for all subjects.
param_initial = cbind(runif(numIter)*param_up[1], runif(numIter)*param_up[2], runif(numIter)*param_up[3])

global_pars = data.frame(rho=NULL, lambda=NULL, tau=NULL, AIC=NULL, BIC=NULL, subjID=NULL)

# for each subject.. 
for (i in 1:N) {
  tmp_ID = allSubjs[i]  # tmpID. e.g., tmpID = 2
  tmp_data = subset(dat, subjID == tmp_ID)  # select only tmp_ID's data
  global_estimate = 1e6   # a very large number so that it can be rejected with an initial MLE estimate
  global_mle = NULL
  
  for (iter in 1:numIter) {
    # fine MLE estimates of the current subject
    mle_ra_prospect <- optim(param_initial[iter, ], ra_prospect_mle, method="L-BFGS-B", 
                             lower=param_low, upper=param_up, tmp_data = tmp_data)
    #print(mle_ra_prospect$value) # for debugging, show outputs
    
    # Replace the results if the latest optimization yields better result
    if(mle_ra_prospect$value < global_estimate) {
      global_estimate <- mle_ra_prospect$value
      global_mle = mle_ra_prospect
    }
  }
  
  global_pars[i, "rho"] = global_mle$par[1]
  global_pars[i, "lambda"] = global_mle$par[2]
  global_pars[i, "tau"] = global_mle$par[3]
  global_pars[i, "AIC"] = 2*global_mle$value + 2*numPars
  global_pars[i, "BIC"] = 2*global_mle$value + numPars*log(T)
  global_pars[i, "subjID"] = tmp_ID
  
  cat("End of modeling subject ID =", tmp_ID, "\n")
}

print(global_pars)

sum_AIC = sum(global_pars$AIC)
print(sum_AIC)

sum_BIC = sum(global_pars$BIC)
print(sum_BIC)


# compare individual MLE and HBA estimates
# first estimate the data with hBayesDM
output_ra_prospect = ra_prospect("example", 2000, 1000, 2, 2)

# compare MLE and HBA
plot(output_ra_prospect$allIndPars$rho, global_pars$rho); abline(0,1)
plot(output_ra_prospect$allIndPars$lambda, global_pars$lambda); abline(0,1)
plot(output_ra_prospect$allIndPars$tau, global_pars$tau); abline(0,1)







