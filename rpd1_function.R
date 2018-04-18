##only for rpd1

#Creating 3*4*128 answer sheet from "correctdeck matrix"
correctdeckmat=read.table("correctdeck matrix.txt",header=1) 
matchstack=array(rep(0,(3*4*128)),dim=c(3,4,128)) 
for (i in 1:4) matchstack[,i,] <- array(as.numeric(correctdeckmat[,]==i),dim=c(3,128))

#Modulating for parameter bounds
stretchpars=function(opars) -log((ub-lb)/(opars-lb)-1)	#opars=original pars
contractpars=function(spars) (ub-lb)/(exp(-spars)+1)+lb		#spars=stretched pars

powerrize=function(pow,tempvec) (tempvec^pow)/sum(tempvec^pow)

#Calculate deck-probabilities(4) based on attention(3)
cattpredp <- matrix(NA, ncol=128, nrow=4)
for(trial in 1:128){
  cattpredp[,trial]=t(c(1/3,1/3,1/3)%*%matchstack[,,trial])
}

cattG2fun=function(temppars,tempchoices){ #attention pars(3)
  templength=length(tempchoices)
  tempchoiceprob=matrix(c(tempchoices==1,tempchoices==2,tempchoices==3,tempchoices==4),nrow=4,ncol=templength,byrow=TRUE)*(cattpredp[,1:length(tempchoices)])
  tempchoiceprob=colSums(tempchoiceprob)
  tempchoiceprob=tempchoiceprob[2:templength]*.9998+.0001		#removes trial 1 and rescales. why?
  return(-2*sum(log(tempchoiceprob)))
}


#Calculate deck-probability(4) under 24 models
#temppars = c(r,p,d,i)



vattpredpfun9=function(temppars,tempchoices,tempreinf){
  r=temppars[1];p=temppars[2];d=temppars[3]
  
  templength=length(tempchoices)
  curatt=rep((1/3),3)
  subjattmat=matrix(nrow=3,ncol=templength)
  predpmat=matrix(-1,nrow=4, ncol=templength)
  subjattmat[,1]=curatt
  for (temptrial in 1:(templength-1))
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