# ------ Preprocess data for friendship network example -------

covData <- read.table("data-raw/comm8_att.dat", skip = 10, stringsAsFactors = FALSE)
names(covData) <- c("Sex", "Race", "Grade")

network <- read.table("data-raw/comm8.dat", skip = 4)

adjMat <- matrix(0, nrow = nrow(covData), ncol = nrow(covData))
for (i in 1:nrow(network)) {
  adjMat[network[i, 1], network[i, 2]] = 1
}
adjMat <- adjMat + t(adjMat) # adjacency matrix undirected
adjMat <- 1 * (adjMat > 1)

#delete missing values
remove <- (covData$Sex == 0) | (covData$Race == 0) | (covData$Grade == 0)
adjMat <- adjMat[!remove, !remove]
covData <- covData [!remove, ]

covData <- data.frame(lapply(covData, as.factor))
levels(covData$Sex) <- c("male","female")
levels(covData$Race) <- c("white","hispanic","mixed/other")

save(adjMat,covData,file = "data/friendshipNet.Rdata")
