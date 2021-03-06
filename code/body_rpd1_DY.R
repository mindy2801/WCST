

rm(list=ls())
source("function_DY.R")

#Set data
datname="managed sample data.txt"
rawdatamat=read.table(datname,encoding="UTF-8")  
subjlabels=rawdatamat[,257]				#reads extra information in datafile if available
subjgroup= rawdatamat[,258]

#Set simulation conditions
maxiter <- 50
maxsubj <- 88
lengthvec <- 128-rowSums(rawdatamat[,1:128]==0)
modelstorun <- 5

parbounds <- c(0,0,.01,.01,1,1,5,5)  #boundaries for r, p, d, i
  lb=parbounds[1:4]
  ub=parbounds[5:8]

stretchpars=function(opars) -log((ub-lb)/(opars-lb)-1)	#opars=original pars
contractpars=function(spars) (ub-lb)/(exp(-spars)+1)+lb		#spars=stretched pars


  
#Generate Model Names

freeparsmat <- expand.grid(r=c("r",""),i=c("i","0","1"),d=c("d",""),p=c("p",""))
  freeparsmat <- as.matrix(freeparsmat)
  freeparsmat <- freeparsmat[,c(1,4,3,2)] #Needed to match the sequence as the original one.

fixedvalsmat <- expand.grid(r=c(-1,1),i=c(-1,0.0001,1-1e-8),d=c(-1,1),p=c(-1,1)) # -1 means free parameter
  fixedvalsmat <- as.matrix(fixedvalsmat)
  fixedvalsmat <- fixedvalsmat[,c(1,4,3,2)]

pequalsrmat <- fixedvalsmat[,"p"]
  pequalsrmat[pequalsrmat==-1] <- 0
modnames <- apply(freeparsmat, 1, paste, collapse="")

#Initailize stacks

twoLLstack <- array(rep(NA,(maxsubj*maxiter*24)),dim=c(maxsubj,maxiter,24)) #row is subj, col is LL, dim is number of models
BICstack=array(rep(NA,(maxsubj*maxiter*24)),dim=c(maxsubj,maxiter,24)) 

parstack <- array(NA, c(maxsubj, 4, 24)) #col is parameters
finalLLstack <- array(NA, dim=c(maxsubj,1,24))
finalBICstack <- array(NA, dim=c(maxsubj,1,24))

##For hypothesized models


#Run simulations

for (cursubj in 50:88){ #for 88 subjects

  
  curlength=lengthvec[cursubj]
  curchoices=data.frame(rawdatamat[cursubj,1:curlength])
  curreinf=data.frame(rawdatamat[cursubj,129:(128+curlength)])
  
  for (curmod in modelstorun){  #for 24 models
    
    temppars <- runif(4)*ub
    setmod <- optim(temppars, vattG2overarchfun, freeletters=freeparsmat[curmod,],fixedvals=fixedvalsmat[curmod,],
                     pequalsr=pequalsrmat[curmod],tempchoices=curchoices,tempreinf=curreinf,predpfun=vattpredpfun9, 
                     method="Nelder-Mead") #abnormal termination happens with L-BFGS-B. have to manually re-range parameters
    
    for (curiter in 1:maxiter){ #run 100 iterations. Optimize the -LL from MLE model

      temppars <- runif(4)*ub
      tempmod <- optim(temppars, vattG2overarchfun, freeletters=freeparsmat[curmod,],fixedvals=fixedvalsmat[curmod,],
                       pequalsr=pequalsrmat[curmod],tempchoices=curchoices,tempreinf=curreinf,predpfun=vattpredpfun9, 
                       method="Nelder-Mead") #abnormal termination happens with L-BFGS-B. have to manually re-range parameters
      
      twoLLstack[cursubj,curiter,curmod] <- tempmod$value
      BICstack[cursubj,curiter,curmod] <- tempmod$value+sum(freeparsmat[curmod,]!="")*log(curlength-1)
        
      if (tempmod$value < setmod$value){ #Stack parameters
        setmod <- tempmod}
      
      roundpars <- round(contractpars(tempmod$par),3)
      print(noquote(c("subj#=",cursubj," iter=",curiter," model=",modnames[curmod], "  -2LL=",round(tempmod$value,3) )))
      print(noquote(c("r=",roundpars[1],"  p=",roundpars[2],"  d=",roundpars[3],"   i=",roundpars[4])))
      print(noquote(""))
      flush.console()
      
      
    } #iteration loop
    
    #Calculate information criteria
    parstack[cursubj,,curmod] <- contractpars(setmod$par)
    finalLLstack[cursubj,,curmod] <- setmod$value
    finalBICstack[cursubj,,curmod] <- tempmod$value+sum(freeparsmat[curmod,]!="")*log(curlength-1)
      
  } #model loop
  
  

} #subject loop
  














##For baseline model
#Calculate information criteria for baseline model. case3로 해보자!
deckbaseG2 <- c()
deckbaseG2_DY <- c()
catt33G2 <- c()
deckbaseBIC <- c()

for(cursubj in 1:maxsubj){
  
  curlength=lengthvec[cursubj]
  curchoices=data.frame(rawdatamat[cursubj,1:curlength])
  curreinf=data.frame(rawdatamat[cursubj,129:(128+curlength)])
  
  deckobsf <- c() 
  for (i in 1:4){deckobsf[i] <- c(sum(curchoices==i))}
  deckexpf <- sum(deckobsf)*c(1/4,1/4,1/4,1/4) #Expected frequency assuming independence
  deckobsp <- deckobsf/lengthvec[cursubj]
  
  deckbaseG2[cursubj]=-2*sum(deckobsf*log(deckobsp)) #original code 지금 이건 loglikelihood랑 g2를 섞은것 같은데...
  deckbaseG2_DY[cursubj] <- -2*lengthvec[cursubj]*log(0.25)
  catt33G2[cursubj]=cattG2fun(rep((1/3),3),curchoices) #G2 아니고 2LL임. attention에 따라 deckchoice probability를 준 뒤, 그 probability를 case3
  deckbaseBIC[cursubj]=deckbaseG2[cursubj]+3*log(curlength-1) #왜 1개 빼지? cattg2fun에서도 그러던데..그럼 deckbaseg2에서도 빼야하는거 아님?
  
  
}


########여기서부터 내 실습
cursubj <- 1

curlength=lengthvec[cursubj]
curchoices=data.frame(rawdatamat[cursubj,1:curlength])
curreinf=data.frame(rawdatamat[cursubj,129:(128+curlength)])

deckbaseG2 <- c()
catt33G2 <- c()
deckbaseBIC <- c()

deckobsf <- c() 
for (i in 1:4){deckobsf[i] <- c(sum(curchoices==i))}
deckexpf <- sum(deckobsf)*c(1/4,1/4,1/4,1/4) #Expected frequency assuming independence
deckobsp <- deckobsf/lengthvec[cursubj]

##Case1. original code
deckbaseG2[cursubj]=-2*sum(deckobsf*log(deckobsp)) #original code 지금 이건 loglikelihood랑 g2를 섞은것 같은데...
catt33G2[cursubj]=cattG2fun(rep((1/3),3),curchoices) #G2 아니고 2LL임. attention에 따라 deckchoice probability를 준 뒤, 그 probability를 case3
deckbaseBIC[cursubj]=deckbaseG2[cursubj]+3*log(curlength-1) #왜 1개 빼지? cattg2fun에서도 그러던데..그럼 deckbaseg2에서도 빼야하는거 아님?


##Case3. choice에 대한 probability sum.  lengthvec*log(0.25) 이게 single trial에 대한 multinomial을 우도함수로 사용한 것인듯?
-2*lengthvec[cursubj]*log(0.25)

##
curmod <- 5
catt33G2 <- array(NA, c(maxsubj, 1, 24)) 
catt33BIC <- array(NA, c(maxsubj, 1, 24)) 
for(cursubj in 1:88){
  curlength=lengthvec[cursubj]
  curchoices=data.frame(rawdatamat[cursubj,1:curlength])
  curreinf=data.frame(rawdatamat[cursubj,129:(128+curlength)])
  
  
  catt33G2[cursubj,1,curmod] <- cattG2fun(rep(1/3,3),curchoices)
  catt33BIC[cursubj,1,curmod] <-catt33G2[cursubj,1,curmod]+3*log(curlength-1)
}


####summarize tables 1:49, 50:88 control/sdi
r_BIC <- finalBICstack[,,5]
r_LL <- finalLLstack[,,5]
r_par <- parstack[,,5]
r_baseLL <- catt33G2[,,5]
r_baseBIc <- catt33BIC[,,5]

mean_par <- rbind(control=colMeans(r_par[1:49,]), sdi=colMeans(r_par[50:88,]))
median_par <- rbind(control=apply(r_par[1:49,],2,median), sdi=apply(r_par[50:88,],2,median))
sd_par <- rbind(control=apply(r_par[1:49,],2,sd), sdi=apply(r_par[50:88,],2,sd))

BIC_5 <- mean(r_BIC)
LL_base <- mean(r_baseLL)
BIC_base <- mean(r_baseBIc)
