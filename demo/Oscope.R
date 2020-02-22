
library(Oscope)

data(OscopeExampleData)
str(OscopeExampleData)
set.seed(100)

#################
# Normalization
#################
Sizes <- MedianNorm(OscopeExampleData)
DataNorm <- GetNormalizedMat(OscopeExampleData, Sizes)

#################
# Select high mean high variance genes
#################
MV <- CalcMV(Data = OscopeExampleData, Sizes = Sizes)
str(MV$GeneToUse)
DataSubset <- DataNorm[MV$GeneToUse,]


#################
# Rescale data
#################
DataInput <- NormForSine(DataNorm)


#################
# Paired sine model
#################
SineRes <- OscopeSine(DataInput)
str(SineRes)


#################
# K-medoids clustering
#################
KMRes <- OscopeKM(SineRes, maxK = 5)
print(KMRes)


#################
# Flag clusters with small within-cluster phase shift
#################
ToRM <- FlagCluster(SineRes,KMRes,DataInput)
print(ToRM$FlagID)
KMResUse <- KMRes[-ToRM$FlagID]

#################
# Extended nearest insertion
#################
ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput, NCThre = 100)
print(ENIRes)


par(mfrow = c(3,2))
for(i in 1:6)
	plot(DataNorm[KMResUse[[1]][i], ENIRes[[1]]],
			 xlab = "Recovered order", ylab = "Expression",
			 main = KMResUse[[1]][i])



par(mfrow = c(3,2))
for(i in 1:6)
	plot(DataNorm[KMResUse[[2]][i], ENIRes[[2]]],
			 xlab = "Recovered order", ylab = "Expression",
			 main = KMResUse[[2]][i])
