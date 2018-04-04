
matchstack=array(rep(0,(3*4*128)),dim=c(3,4,128)) 
correctdeckmat=read.table("correctdeck matrix.txt",header=1)
for (i in 1:4) matchstack[,i,]=array(as.numeric(correctdeckmat[,]==i),dim=c(3,128))

stretchpars=function(opars) -log((ub-lb)/(opars-lb)-1)	#opars=original pars
contractpars=function(spars) (ub-lb)/(exp(-spars)+1)+lb		#spars=stretched pars

scale3to2=function(temppars) c(temppars[1],temppars[2]/(1-temppars[1]))
#this takes a vector of 3 p's that would sum to 1 and rescales them into 2 pars that range from 0-1
#this is done so optimization with range 0-1 can work.  e.g. (.6,.1,.3) -> (.6,.25)

scale2to3=function(temppars)  c(temppars[1],(1-temppars[1])*temppars[2],1-temppars[1]+(temppars[1]-1)*temppars[2])
#this takes a vector of 2 pars that range from 0-1 and transforms back to 3 p's that sum to 1

powerrize=function(pow,tempvec) (tempvec^pow)/sum(tempvec^pow)
#takes a vector, raises to a power, then rescales to sum to 1

cattpredpfun=function(temppars,templength) 
  #predicted probabilities assuming Constant Attention weights
  #temppars is a 3 parameter pars vector representing the attention weights to color, form, and number
{
  tempmat=matrix(-1,ncol=templength,nrow=4)
  for (temptrial in 1:templength) tempmat[,temptrial]=t(temppars%*%matchstack[,,temptrial])
  tempmat
}

#vattpredpfun9 is an overarching all purpose model function calculating the predicted probability
#takes parameters r, p, d, i ("i" might be called "g" or "f" in paper)
#when trial>m the "aha" experience has occurred. globalm is a global variable
#freeletters is a vector of characters in any order for the free parameter letters
#fixedvals is an ordered vector numbers that only matters for nonfree parameters
#pequalsr is boolean for whether constraint p=r is true
vattpredpfun9=function(temppars,freeletters,fixedvals,pequalsr,tempchoices,tempreinf) 
{
  if(sum(freeletters=="r")>0) r=temppars[1] else r=fixedvals[1]
  if(sum(freeletters=="p")>0) p=temppars[2] else p=fixedvals[2] 
  if(sum(freeletters=="d")>0) d=temppars[3] else d=fixedvals[3]
  if(sum(freeletters=="i")>0) i=temppars[4] else i=fixedvals[4]
  if(pequalsr) p=r
  
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
  return(predpmat)
}

cattG2fun=function(pars2,tempchoices) #generates G2 for constant attention
  #this takes the 2 parameter pars vector (because only two are free)
{
  temppars=scale2to3(pars2)
  templength=length(tempchoices)
  tempchoiceprob=matrix(c(tempchoices==1,tempchoices==2,tempchoices==3,tempchoices==4),nrow=4,ncol=templength,byrow=TRUE)*(cattpredpfun(temppars,templength))
  tempchoiceprob=colSums(tempchoiceprob)
  tempchoiceprob=tempchoiceprob[2:templength]*.9998+.0001		#removes trial 1 and rescales
  return(-2*sum(log(tempchoiceprob)))
}

vattG2overarchfun=function(temppars,tempparbounds,freeletters,fixedvals,pequalsr,tempchoices,tempreinf,predpfun) 
{
  origpars=contractpars(temppars)
  templength=length(tempchoices)
  tempchoiceprob=matrix(c(tempchoices==1,tempchoices==2,tempchoices==3,tempchoices==4),nrow=4,ncol=templength,byrow=TRUE)*( predpfun(origpars,freeletters,fixedvals,pequalsr,tempchoices,tempreinf))
  tempchoiceprob=colSums(tempchoiceprob)
  tempchoiceprob=tempchoiceprob[2:templength]*.9998+.0001		#removes trial 1 and rescales
  return(-2*sum(log(tempchoiceprob)))
}