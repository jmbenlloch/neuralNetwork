setwd("/home/jmbenlloch/next/projects/petalo/neuralNetworks")
data <- as.matrix(read.table("logNeuralNet.txt"))

data[,602:604]

data[,602] <- (data[,602]+50)/100

n<- nrow(data)
nent<- ceiling(0.7*n)
ntest<- n-nent
indin<- 1:n
indient<- sort(sample(indin,nent))
inditest<- setdiff(indin,indient)
dataEnt <- data[indient,] ; dataTest <- data[inditest,]

library(nnet)
#red<- nnet(x=dataEnt[,2:601],y=dataEnt[,602],size=100,maxit=5000,entropy=T,MaxNWts=100000)
red<- nnet(x=dataEnt[,2:101],y=dataEnt[,602],size=100,maxit=5000,entropy=T,MaxNWts=100000)

red
summary(red)
predict(red,dataTest[2,2:601])
testResults <- (predict(red,dataTest[,2:601]) - dataTest[,602])*100
trainResults <- (predict(red,dataEnt[,2:601]) - dataEnt[,602])*100
hist(testResults,breaks=50)
hist(traintResults,breaks=50)

write.table(testResults, "testResults.txt", sep="\t", row.names=FALSE) 
write.table(trainResults, "trainResults.txt", sep="\t", row.names=FALSE) 

