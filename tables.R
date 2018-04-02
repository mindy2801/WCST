##Table.2 mean BIC for 24 models
meanBIC <- colMeans(BICfinalmat); meanBIC
r=rep(c("free","fixed"),12); r <- as.factor(r)
p=c(rep("free",12),rep("fixed",12)); p <- as.factor(p)
d=rep(c(rep("free",6),rep("fixed",6)),2); d <- as.factor(d)
f=rep(c(rep("free",2),rep(0,2),rep(1,2)),4); f <- as.factor(f)
t2 <- data.frame(meanBIC, r, p, d, f)
t2_nr <- t2[t2$r=="free",]; t2_nr <- t2_nr[,-2]
vars <- t2_nr[,1]
v1 <- c(vars[c(1,4,7,10)])
v2 <- c(vars[c(2,5,8,11)])
v3 <- c(vars[c(3,6,9,12)])
vars <- data.frame(v1, v2, v3)
colnames(vars)<-c("f free", "f = 0", "f = 1")

# https://www.r-bloggers.com/fast-track-publishing-using-knitr-table-mania-part-iv/
library(xtable)
alt_vars <- cbind(c("p free", "", "p = r", ""), 
                  c("d free", "d = 1", "d free", "d = 1")
                  )
colnames(alt_vars) <- gsub("", "\n", colnames(alt_vars))
options(xtable.html.table.attributes = 
          list(style=sprintf("style='%s'",
                             paste("border:0",
                                   "border-top: 1px solid grey", 
                                   "border-bottom: 1px solid grey",
                                   sep="; "))))
print(xtable(alt_vars, caption="An xtable example"), type="html", include.rownames = FALSE)

##Fig.2 observed scores for categories completed/presev error/nonperserv error/set failures


##Fig.2 model simulation

##Table.3.1 mean parameters for rpd1 model
meanPar <- colSums(parstack[,1:3,5]); meanPar <- t(meanPar)
library(knitr)
kable(meanPar)
##Table.3.2
control <- parstack[1,1:3,5]; control
SDI <- colMeans(parstack[2:3,1:3,5]); SDI
kable(t3.2)
###
col <- sprintf("Cell %d:%%d", 1:9)
vars <- sapply(1:3, function(i) sprintf(col, i))
alt_vars <- cbind(
  Rowgrp1 = c("Mjr group 1", "", "", 
              "Mjr group 2", "", "", "", "", ""),
  Rowgrp2 = c("Group 1", "", "", 
              "Group 2", "", 
              "Group 3", "", "", ""),
  Rownames= rownames(vars), 
  vars)
colnames(alt_vars) <- gsub("
                           ", "\n", colnames(alt_vars))
# rownames(vars) <- NULL
options(xtable.html.table.attributes = 
          list(style=sprintf("style='%s'",
                             paste("border:0",
                                   "border-top: 1px solid grey", 
                                   "border-bottom: 1px solid grey",
                                   sep="; "))))
print(xtable(alt_vars, caption="An xtable example"), type="html", include.rownames = FALSE)