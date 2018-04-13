#For hypothesized model

  #For 88 subjects

    #100 iterations
      #MLE model that returns sum of -LL
      #Optimize the -LL from MLE model

    #Stack parameters
    #Calculate BIC


########

#For baseline model
rm(list=ls())
source("probability.R")

#Set data
datname="WCST Sample data.txt"
rawdatamat=read.table(datname,encoding="UTF-8")  
subjlabels=rawdatamat[,257]				#reads extra information in datafile if available
subjgroup= rawdatamat[,258]

#Set simulation conditions
maxiter <- 10
maxsubj <- nrow(rawdatamat)
lengthvec <- 128-rowSums(rawdatamat[,1:128]==0)

parbounds <- c(0,0,.01,.01,1,1,5,5)  #boundaries for r, p, d, i
  lb=parbounds[1:4]
  ub=parbounds[5:8]

#Generate Model Names

freeparsmat <- expand.grid(r=c("r",""),i=c("i","0","1"),d=c("d",""),p=c("p",""))
  freeparsmat <- freeparsmat[,c(1,4,3,2)] #Needed to match the sequence as the original one.

fixedvalsmat <- expand.grid(r=c(-1,1),i=c(-1,0.0001,1-1e-8),d=c(-1,1),p=c(-1,1)) # -1 means free parameter
  fixedvalsmat <- fixedvalsmat[,c(1,4,3,2)]

pequalsrmat <- fixedvalsmat$p
  modnames <- apply(freeparsmat, 1, paste, collapse="")

#Initailize stacks

twoLLstack <- array(NA, c(maxsubj, 1, 1)) #row is subj, col is LL, dim is number of models
parstack <- array(NA, c(maxsubj, 4, 1)) #col is parameters


#Calculate information criteria for baseline model

curlength=lengthvec[cursubj]
curchoices=data.frame(rawdatamat[cursubj,1:curlength])
curreinf=data.frame(rawdatamat[cursubj,129:(128+curlength)])

deckobsf <- c() 
  for (i in 1:4){deckobsf[i] <- c(sum(curchoices==i))}
deckexpf <- sum(deckobsf)*c(1/4,1/4,1/4,1/4) #Expected frequency based on independence
2*sum(deckobsf*log(deckobsf/deckexpf))

deckobsp=c(mean(curchoices==1),mean(curchoices==2),mean(curchoices==3),mean(curchoices==4))
deckobsf=deckobsp*curlength
expected.count <- sum(deckobsf)*c(1/4,1/4,1/4,1/4)

##NEED TO CHECK 멀티노미얼 이용해서 likelihood 직접 계산해야하는 거 아닌가?????????????
deckbaseG2[cursubj]=-2*sum(deckobsf*log(deckobsp)) #original code 지금 이건 loglikelihood랑 g2를 섞은것 같은데...

deckobsf <- c()  #What I think
for (i in 1:4){deckobsf[i] <- c(sum(curchoices==i))}
deckexpf <- sum(deckobsf)*c(1/4,1/4,1/4,1/4) #Expected frequency based on independence 
2*sum(deckobsf*log(deckobsf/deckexpf)) #이게 모든 덱을 랜덤하게 선택했을 때를 영가설로 잡은 g2 아닌가?
sum(deckobsf/lengthvec[cursubj])

##

catt33G2[cursubj]=cattG2fun(scale3to2(rep((1/3),3)),curchoices)
deckbaseBIC[cursubj]=deckbaseG2[cursubj]+3*log(curlength-1)


#Run simulations
for (subj in 1:maxsubj){

  temppars <- runif(4)*ub
  
  curlength=lengthvec[subj]
  curmod <- 5 #
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