dataPath = system.file("extdata/prl_exampleData.txt", package="hBayesDM")
exmp <- read.table(dataPath, header=T)
str(exmp); head(exmp)
  #20 Subjects, 1-100 Trials, 1-2 choice, 1 or -1 outcome
