
#This program performs sequential learning models for the Wisconsin Card Sort Task.
#Programmed by Anthony Bishara. 
#More details can be found in:
#Bishara, A. J., Kruschke, J. K., Stout, J. C., Bechara, A., McCabe, D. P., & Busemeyer, J. R. (2010). 
#Sequential learning models for the Wisconsin Card Sort Task...Journal of Mathematical Psychology.

#Code 17b, June 2012 (clearer instructions/save names, more robust handling of errors)


#Instructions:
#1) Run R. 
#2) In the program, click File, then Change Dir., then choose the directory your data is in.
#3) The data file should be in a tab-delimited text file (UTF-8 encoding) with 1 row per subject. 
#The first 128 columns should list deck choice, with 1=far left deck, 4=far right, and 0 indicating empty
#trials (e.g. if all 6 categories are completed by trial 110, columns 111-128 would consist of zeros).  
#Columns 129-256 should show 1 for correct and 0 for incorrect. Values of -1 are ignored.  
#You can use columns 257 and beyond to put any other information in you want (e.g., the subject number).  
#THERE SHOULDN'T BE ANY ROW OR COLUMN LABELS.  
#The current program was built for standard WCST administration with up to 128 trials, 
#the 10-in-a-row criteria for category achievement, 6 categories to be achieved, and the 
#standard ordering of cards presented.
#4) Download or create a file to represent the correct answers in the WCST version that you used.
#The web page where you found this code contains a file called "correctdeck matrix.txt".  
#That file is for the standard, Heaton version with the standard card ordering.
#This file should placed in the same directory as your code and data.
#5) Modify the  section of the code below called "### MODIFY THIS SECTION ###"
#6) Select all the text in this document (ctrl+a), then copy and paste into R
#7) The output file is labelled "Code17b_" plus details about the modeling.   


rm(list=ls(all=TRUE))  					#clears memory

### MODIFY THIS SECTION ###
datname="WCST Sample data.txt"
maxiter=100 				#Number of starting parameters to iterate through; default is 100
modelstorun=5    #lists the nested models that'll be run; #5 is the default, best fitting model in article; use "modelstorun=1:24" to run all (slower)
parbounds=c(0,0,.01,.01,1,1,5,5)  #lower boundaries for parameters r, p, d, i and upper boundaries for parameters r, p, d, i
#############################


if(sum(rownames(installed.packages())=="fOptions")==0)   #If the fOptions package is not installed, .
  install.packages("fOptions") 				#.install it now.  
library(fOptions)  						#load fOptions into memory.

rawdatamat=read.table(datname,encoding="UTF-8")  
#If you receive an error after this line, then the datafile isn't being read in properly.
#One possible problem is that the datafile doesn't have UTF-8 Encoding.  
#One solution is to use the Notepad program to resave the datafile.
#At the bottom of the Save dialog window, next to "Encoding:", choose "UTF-8" before saving the datafile.
#If that does not fix the problem, try opening the file in MS Word (using windows default encoding), 
#delete the unusual characters, and then resave it (again using windows default encoding).

subjlabels=rawdatamat[,257]				#reads extra information in datafile if available
subjgroup= rawdatamat[,258]

#higher values will lead to more accurate but slower processing


lb=parbounds[1:4]			#global variables for bounds (lb=lower bounds, ub=upper bounds)
ub=parbounds[5:8]

itercolumns=maxiter
subjectsmodeled=1:length(rawdatamat[,1]) #subject number
savelabel=paste("Code17b_",strtrim(datname,6),"_s",min(subjectsmodeled),"-",max(subjectsmodeled),"_iter",maxiter,sep="") #saved filename

#define a set of 128 match matrices (3 dimensions x 4 decks x 128 trials)
matchstack=array(rep(0,(3*4*128)),dim=c(3,4,128)) 
correctdeckmat=read.table("correctdeck matrix.txt",header=1) #Using matrix.txt, creating anwser sheet
for (i in 1:4) matchstack[,i,]=array(as.numeric(correctdeckmat[,]==i),dim=c(3,128)) #IMPORTANT TO UNDERSTAND

stretchpars=function(opars) -log((ub-lb)/(opars-lb)-1)	#opars=original pars #WHAT IS IT DOING?
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

#PROGRAM BODY
numsubj=nrow(rawdatamat)
deckbaseG2=rep(-1,numsubj)	#deckbase is a baseline model which predicts .25 probability of choosing each deck
deckbaseBIC=rep(-1,numsubj)
catt33G2=rep(-1,numsubj) #catt33 is a different baseline model which has constant .33 attention to each dimension
lengthvec=128-rowSums(rawdatamat[,1:128]==0) #number of trials

freeparsmat=matrix(nrow=24,ncol=4,dimnames=list(NULL,c("r","p","d","i")))  #setting up a matrix that lists the free pars of each of 24 models
fixedvalsmat=matrix(-1,nrow=24,ncol=4,dimnames=list(NULL,c("r","p","d","i"))) #parameter constraint values of 24 models (when -1, free parameter)
pequalsrmat=matrix(0,nrow=24,ncol=1)
#until now, empty vectors are created

for (ploop in 0:1) for (dloop in 0:1) for (iloop in 0:2) for (rloop in 0:1)
{
  rowloop=ploop*12 + dloop*6 + iloop*2 + rloop +1
  freeparsmat[rowloop,1]=(if(rloop==0) "r" else "")
  freeparsmat[rowloop,2]=(if(ploop==0) "p" else "")
  freeparsmat[rowloop,3]=(if(dloop==0) "d" else "")
  freeparsmat[rowloop,4]=(if(iloop==0) "i" else "")
  if(dloop==1) fixedvalsmat[rowloop,3]=1-1e-8 #WHAT IS THIS NUMBER?
  if(iloop==1) fixedvalsmat[rowloop,4]=.0001
  if(iloop==2) fixedvalsmat[rowloop,4]=1
  if(rloop==1) fixedvalsmat[rowloop,1]=1
  if(ploop==1) pequalsrmat[rowloop,1]=1 #WHY NOT fixedvalsmat?
}

parnames=apply(freeparsmat,1,paste,collapse="")
modnames=parnames
modnames[freeparsmat[,"i"]==""]=paste(modnames[freeparsmat[,"i"]==""], round(fixedvalsmat[freeparsmat[,"i"]=="","i"],1),sep="")
  #when "i" colmun(4th) is empty, at that location within modnames, 1-digit round of feixedvalsmat's 4th column is pasted

G2stack=array(rep(NA,(numsubj*itercolumns*24)),dim=c(numsubj,itercolumns,24)) #3*100*24, by each iteration(till 100) maybe store estimated parameter for each subject with 24 models?
BICstack=array(rep(NA,(numsubj*itercolumns*24)),dim=c(numsubj,itercolumns,24)) 
dimnames(G2stack)=list(paste("s",subjlabels,sep=""),paste("i",1:itercolumns,sep=""),paste("m",1:24,"_",modnames,sep=""))
dimnames(BICstack)=list(paste("s", subjlabels,sep=""),paste("i",1:itercolumns,sep=""),paste("m",1:24,"_",modnames,sep=""))

parstack=array(rep(NA,(numsubj*4*24)),dim=c(numsubj,4,24)) #3*4*24
dimnames(parstack)=list(paste("s",subjlabels,sep=""),c("r","p","d","i"),paste("m",1:24,"_",modnames,sep=""))
presobelmat=sobelmat=runif.sobol(maxiter,4,1) #Uniform scrambled sobol sequence. WHAT IS SOBOL SEQUENCE?
for (i in 1:4) sobelmat[,i]=presobelmat[,i]*(parbounds[i+4]-parbounds[i])+parbounds[i] #Using sobol number, recreating numbers within bound

starttime=Sys.time()
for (cursubj in subjectsmodeled)	#loop across subjects
{
  curlength=lengthvec[cursubj]
  curchoices=data.frame(rawdatamat[cursubj,1:curlength]) #put out valid choices
  curreinf=data.frame(rawdatamat[cursubj,129:(128+curlength)]) #put out correct or not
  deckobsp=c(mean(curchoices==1),mean(curchoices==2),mean(curchoices==3),mean(curchoices==4)) #mean of each choices
  deckobsf=deckobsp*curlength #frequency of each choices
  deckbaseG2[cursubj]=-2*sum(deckobsf*log(deckobsp)) #WHAT IS IT DOING?
  catt33G2[cursubj]=cattG2fun(scale3to2(rep((1/3),3)),curchoices) #WHAT IS IT DOING???? Why remove trial 1 and rescales?
  
  for (curmod in modelstorun)	#loop across models (different parameter constraints)
  {
    for (curiter in 1:itercolumns)
    curiter=0
    contiter=TRUE
    while(contiter)		#loop across different iterations (i.e., starting parameters)
    {
      curiter=curiter+1
      pars4init=sobelmat[curiter,]
      spars4init=stretchpars(pars4init) #Why is stretching needed?
      tempmod=optim(spars4init,vattG2overarchfun,tempparbounds=parbounds,freeletters=freeparsmat[curmod,],fixedvals=fixedvalsmat[curmod,],pequalsr=pequalsrmat[curmod,],tempchoices=curchoices,tempreinf=curreinf,predpfun=vattpredpfun9,method="Nelder-Mead")
      G2stack[cursubj,curiter,curmod]=tempmod$value #WHAT IS OPTIM FUNCTION?
      roundpars=round(contractpars(tempmod$par),3)					
      print(noquote(c("subj#=",cursubj," iter=",curiter," model=",modnames[curmod], "  -2LL=",round(tempmod$value,3) )))
      print(noquote(c("r=",roundpars[1],"  p=",roundpars[2],"  d=",roundpars[3],"   i=",roundpars[4])))
      print(noquote(""))
      flush.console()
      
      if(curiter==1) parstack[cursubj,,curmod]=contractpars(tempmod$par) else
      {
        if(tempmod$value<min(G2stack[cursubj,1:curiter-1,curmod])) parstack[cursubj,,curmod]=contractpars(tempmod$par)
      }
      BICstack[cursubj,curiter,curmod]=G2stack[cursubj,curiter,curmod]+sum(freeparsmat[curmod,]!="")*log(curlength-1)
      if(curiter>=maxiter) contiter=FALSE
    }
  }
  deckbaseBIC[cursubj]=deckbaseG2[cursubj]+3*log(curlength-1)
  
  
}
timeelapsed=Sys.time()-starttime
timeelapsed

G2finalmat=array(rep(NA,(length(subjectsmodeled)*24)),dim=c(length(subjectsmodeled),24))
BICfinalmat=array(rep(NA,(length(subjectsmodeled)*24)),dim=c(length(subjectsmodeled),24))
dimnames(G2finalmat)=list(paste("s",subjlabels[subjectsmodeled],sep=""),paste("G2m",1:24,"_",modnames,sep="")) 
dimnames(BICfinalmat)=list(paste("s",subjlabels[subjectsmodeled],sep=""),paste("BICm",1:24,"_",modnames,sep="")) 
parfinalmat=subjlabels[subjectsmodeled]
parlabels="subjnum"
for (i in modelstorun)
{
  G2finalmat[,i]=apply((G2stack[subjectsmodeled,,i]),1,min,na.rm=1)
  BICfinalmat[,i]=apply((BICstack[subjectsmodeled,,i]),1,min,na.rm=1)
  parfinalmat=cbind(parfinalmat,parstack[subjectsmodeled,,i])
  parlabels=c(parlabels,paste("m",i,"_",modnames[i],"_",freeparsmat[1,],sep=""))
}
G2BICmat=cbind(subjlabels[subjectsmodeled],catt33G2[subjectsmodeled],G2finalmat,BICfinalmat)
colnames(G2BICmat)[1:2]=c("subjnum","catt33G2")
colnames(parfinalmat)=parlabels
write.table(cbind(G2BICmat,parfinalmat), paste(savelabel,".txt",sep=""))
save.image(paste(savelabel,".Rdata",sep=""))
