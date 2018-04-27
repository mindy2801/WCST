

data <- read.table("managed sample data.txt")
library(dplyr)
data <- rename(data, subjID=V257, group=V258)

nsubj <- nrow(data)

#Seperate deck choice and outcome
empty <- matrix(NA, nrow=2*nsubj, ncol=128)

for(j in 1:nsubj){
  empty[2*j-1,1:128] <- as.numeric(data[j,1:128])
  empty[2*j,1:128] <- as.numeric(data[j,129:256])
  }

#Transpose and row-bind trials
t <- c()
for(i in 1:nsubj) {
  t <- rbind(t,  t(empty[(2*i-1):(2*i),]))
}

#Column bind indices
subjID <- rep(data$subjID, each=128)
group <- rep(data$group, each=128)
trial <- rep(1:128, nsubj)
data <- data.frame(t,subjID,group, trial)  

#Remove empty trials
loc <- which(data[,1]==0)
data <- data[-loc,]

data <- data %>% rename(deck=X1, outcome=X2) %>% select(subjID, group, trial, deck, outcome) %>% arrange(subjID)
write.table(data, file="reorganized data.txt")


