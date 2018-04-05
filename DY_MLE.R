#For hypothesized model

  #For 88 subjects

    #100 iterations
      #MLE model that returns sum of -LL
      #Optimize the -LL from MLE model

    #Stack parameters
    #Calculate BIC



#For baseline model
rm(list=ls())
source("probability.R")

datname="WCST Sample data.txt"
rawdatamat=read.table(datname,encoding="UTF-8")  
subjlabels=rawdatamat[,257]				#reads extra information in datafile if available
subjgroup= rawdatamat[,258]

maxiter <- 10
maxsubj <- nrow(rawdatamat)

parbounds=c(0,0,.01,.01,1,1,5,5)  #boundaries for r, p, d, i
  lb=parbounds[1:4]
  ub=parbounds[5:8]

freeparsmat=matrix(nrow=24,ncol=4,dimnames=list(NULL,c("r","p","d","i")))  #setting up a matrix that lists the free pars of each of 24 models
fixedvalsmat=matrix(-1,nrow=24,ncol=4,dimnames=list(NULL,c("r","p","d","i"))) #parameter constraint values of 24 models
pequalsrmat=matrix(0,nrow=24,ncol=1)
  for (ploop in 0:1) for (dloop in 0:1) for (iloop in 0:2) for (rloop in 0:1)
  {
    rowloop=ploop*12 + dloop*6 + iloop*2 + rloop +1
    freeparsmat[rowloop,1]=(if(rloop==0) "r" else "")
    freeparsmat[rowloop,2]=(if(ploop==0) "p" else "")
    freeparsmat[rowloop,3]=(if(dloop==0) "d" else "")
    freeparsmat[rowloop,4]=(if(iloop==0) "i" else "")
    if(dloop==1) fixedvalsmat[rowloop,3]=1-1e-8
    if(iloop==1) fixedvalsmat[rowloop,4]=.0001
    if(iloop==2) fixedvalsmat[rowloop,4]=1
    if(rloop==1) fixedvalsmat[rowloop,1]=1
    if(ploop==1) pequalsrmat[rowloop,1]=1
  }
  
  parnames=apply(freeparsmat,1,paste,collapse="")
  modnames=parnames
  modnames[freeparsmat[,"i"]==""]=paste(modnames[freeparsmat[,"i"]==""], round(fixedvalsmat[freeparsmat[,"i"]=="","i"],1),sep="")
  
  
twoLLstack <- array(NA, c(maxsubj, 1, 1))
parstack <- array(NA, c(maxsubj, 4, 1))
lengthvec=128-rowSums(rawdatamat[,1:128]==0)

for (subj in 1:maxsubj){

  temppars <- runif(4)*ub
  
  curlength=lengthvec[subj]
  curmod <- 5
  curchoices=data.frame(rawdatamat[subj,1:curlength])
  curreinf=data.frame(rawdatamat[subj,129:(128+curlength)])
  
  temppars <- runif(4)*ub
  setmod <- optim(temppars, vattG2overarchfun, freeletters=freeparsmat[curmod,],fixedvals=fixedvalsmat[curmod,],pequalsr=pequalsrmat[curmod,],
                  tempchoices=curchoices,tempreinf=curreinf,predpfun=vattpredpfun9, lower=lb, upper=ub, method="L-BFGS-B")
  
  
    for (iter in 1:maxiter){ 
      temppars <- runif(4)*ub
      tempmod <- optim(temppars, vattG2overarchfun, freeletters=freeparsmat[curmod,],fixedvals=fixedvalsmat[curmod,],pequalsr=pequalsrmat[curmod,],
                       tempchoices=curchoices,tempreinf=curreinf,predpfun=vattpredpfun9, lower=lb, upper=ub, method="L-BFGS-B")
      
      if (tempmod$value < setmod$value){
        setmod <- tempmod
      }
      twoLLstack[subj,,] <- setmod$value
      parstack[subj,,] <- setmod$par
    }
    
}
  

optim(temppars, vattG2overarchfun, freeletters=freeparsmat[curmod,],fixedvals=fixedvalsmat[curmod,],pequalsr=pequalsrmat[curmod,],
                  tempchoices=curchoices,tempreinf=curreinf,predpfun=vattpredpfun9,lower=lb, upper=ub, method="L-BFGS-B") 

optim(temppars, vattG2overarchfun, freeletters=freeparsmat[curmod,],fixedvals=fixedvalsmat[curmod,],pequalsr=pequalsrmat[curmod,],
      tempchoices=curchoices,tempreinf=curreinf,predpfun=vattpredpfun9, method="Nelder-Mead") 

k <-
N <- #subj number?
tBIC <- colSums(twoLLstack) + log(N)*k
tpar <- colMeans(parstack)