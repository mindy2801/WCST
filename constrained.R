#only for rpd1

rm(list=ls())
source("rpd1_function.R")

#Set data
datname="managed sample data.txt"
rawdatamat=read.table(datname,encoding="UTF-8")  
subjlabels=rawdatamat[,257]				#reads extra information in datafile if available
subjgroup= rawdatamat[,258]

#Set simulation conditions
maxiter <- 5
maxsubj <- nrow(rawdatamat)
lengthvec <- 128-rowSums(rawdatamat[,1:128]==0)
modelstorun <- 5

parbounds <- c(0,0,.01,1,1,5)  #boundaries for r, p, d
lb=parbounds[1:3]
ub=parbounds[4:6]

stretchpars=function(opars) -log((ub-lb)/(opars-lb)-1)	#opars=original pars
contractpars=function(spars) (ub-lb)/(exp(-spars)+1)+lb		#spars=stretched pars



#Generate Model Names


#Initailize stacks

twoLLstack <- array(rep(NA,(maxsubj*maxiter*1)),dim=c(maxsubj,maxiter,1)) #row is subj, col is LL, dim is number of models
BICstack=array(rep(NA,(maxsubj*maxiter*1)),dim=c(maxsubj,maxiter,1)) 

parstack <- array(NA, c(maxsubj, 3, 1)) #col is parameters
finalLLstack <- array(NA, dim=c(maxsubj,1,1))
finalBICstack <- array(NA, dim=c(maxsubj,1,1))

##For hypothesized models



for (cursubj in 1:maxsubj){ #for 88 subjects
  
  
  curlength=lengthvec[cursubj]
  curchoices=data.frame(rawdatamat[cursubj,1:curlength])
  curreinf=data.frame(rawdatamat[cursubj,129:(128+curlength)])
  
  
    
    temppars <- runif(3)*ub
    setmod <- optim(temppars, vattG2overarchfun, tempchoices=curchoices,tempreinf=curreinf,predpfun=vattpredpfun9, 
                    method="Nelder-Mead") #abnormal termination happens with L-BFGS-B. have to manually re-range parameters
    
    for (curiter in 1:maxiter){ #run 100 iterations. Optimize the -LL from MLE model
      
      temppars <- runif(3)*ub
      tempmod <- optim(temppars, vattG2overarchfun, tempchoices=curchoices,tempreinf=curreinf,predpfun=vattpredpfun9, 
                       method="Nelder-Mead") #abnormal termination happens with L-BFGS-B. have to manually re-range parameters
      
      twoLLstack[cursubj,curiter,1] <- tempmod$value
      BICstack[cursubj,curiter,1] <- tempmod$value+3*log(curlength-1)
      
      if (tempmod$value < setmod$value){ #Stack parameters
        setmod <- tempmod}
      
      roundpars <- round(contractpars(tempmod$par),3)
      print(noquote(c("subj#=",cursubj," iter=",curiter,"  -2LL=",round(tempmod$value,3) )))
      print(noquote(c("r=",roundpars[1], "  p=",roundpars[2],"  d=",roundpars[3])))
      print(noquote(""))
      flush.console()
      
      
    } #iteration loop
    
    #Calculate information criteria
    parstack[cursubj,,1] <- contractpars(setmod$par)
    finalLLstack[cursubj,,1] <- setmod$value
    finalBICstack[cursubj,,1] <- tempmod$value+3*log(curlength-1)
    
  
  
  
} #subject loop

catt33G2 <- array(NA, c(maxsubj, 1, 1)) 
catt33BIC <- array(NA, c(maxsubj, 1, 1)) 
for(cursubj in 1:88){
  curlength=lengthvec[cursubj]
  curchoices=data.frame(rawdatamat[cursubj,1:curlength])
  curreinf=data.frame(rawdatamat[cursubj,129:(128+curlength)])
  
  
  catt33G2[cursubj,1,1] <- cattG2fun(rep(1/3,3),curchoices)
  catt33BIC[cursubj,1,1] <-catt33G2[cursubj,1,1]+3*log(curlength-1)
}


##

mean_par <- rbind(control=colMeans(parstack[1:49,,]), sdi=colMeans(parstack[50:88,,]))
median_par <- rbind(control=apply(parstack[1:49,,],2,median), sdi=apply(parstack[50:88,,],2,median))
sd_par <- rbind(control=apply(parstack[1:49,,],2,sd), sdi=apply(parstack[50:88,,],2,sd))

BIC_rpd1 <- mean(finalBICstack)
LL_base <- mean(catt33G2)
BIC_base <- mean(catt33BIC)
