#For hypothesized model

  #For 88 subjects

    #100 iterations
      #MLE model that returns sum of -LL
      #Optimize the -LL from MLE model

    #Stack parameters
    #Calculate BIC



#For baseline model
rm(list=ls(all=TRUE))
source("probabiliry.R")

datname="choices correct Bechara SDIandcontrols.txt"
rawdatamat=read.table(datname,encoding="UTF-8")  
subjlabels=rawdatamat[,257]				#reads extra information in datafile if available
subjgroup= rawdatamat[,258]

maxiter <- 100
maxsubj <- nrow[rawdatamat]

parbounds=c(0,0,.01,.01,1,1,5,5)  #boundaries for r, p, d, i
  lb=parbounds[1:4]
  ub=parbounds[5:8]


  
  
  for (subj in 1:maxsubj){

  
  twoLLstack <- array(NA, c(maxsubj, 1, 1))
  parstack <- array(NA, c(maxsubj, 4, 1))
  setmod <- optim(temppars, MLEFUNTION, ... )
  
    for (iter in 1:maxiter){ 
      temppars <- runif(4)*ub
      tempmod <- optim(temppars, MLEFUNCTION, ... )
      
      if (curmod$value < setmod$value){
        setmod <- tempmod
      }
      twoLLstack[subj,,] <- setmod$value
      parstack[subj,,] <- setmod$par
    }
    
    
  }
  
k <-
N <- #subj number?
tBIC <- colSums(twoLLstack) + log(N)*k
tpar <- colMeans(parstack)