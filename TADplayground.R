# this script just plays around with various parameters to make the node graph more recognizable
# and ideally find some parameters that distinguish node graphs between ESC and Neuron


# we'll need these packages later
install.packages("devtools", repos="http://cran.rstudio.com/")
library(devtools)
devtools::install_github("paultpearson/TDAmapper")
library(TDAmapper)
install.packages("fastcluster", repos="http://cran.rstudio.com/")
require(fastcluster)
install.packages("igraph", repos="http://cran.rstudio.com/")
library(igraph)
library(ggplot2)

rawESCdata <- read.csv("~/Dropbox/DS/wTrush/rawESCdata.csv", header=TRUE, stringsAsFactors = FALSE)
rawNeurondata <- read.csv("~/Dropbox/DS/wTrush/rawNeurondata.csv", header=TRUE, stringsAsFactors = FALSE)

########################################################        
##############    chr1 ESC/Neuron            ###########
########################################################        
ESCchr1col <- which(names(rawESCdata)=="chr1")
ESCchr1 <- rawESCdata[2:length(rawESCdata[,ESCchr1col]), (ESCchr1col+1):(ESCchr1col+3)]

ESCchr1[,1] <- as.numeric(ESCchr1[,1])
ESCchr1[,2] <- as.numeric(ESCchr1[,2])
ESCchr1[,3] <- as.numeric(ESCchr1[,3])

names(ESCchr1) = c("startTAD", "endTAD", "interaction")

ESCchr1 <- ESCchr1[1:sum(!is.na(ESCchr1$startTAD)), ]

ESCchr1$dist = 1/ESCchr1$interaction

# now we want to create a distance metric matrix
n1 = max(ESCchr1$endTAD, ESCchr1$startTAD) - (min(ESCchr1$endTAD, ESCchr1$startTAD) - 1)
ESCchr1metric <- matrix(, n1, n1)

for (i in 1:length(ESCchr1$startTAD)) {
        ESCchr1metric[ESCchr1$startTAD[i], ESCchr1$endTAD[i]] = ESCchr1$dist[i]
}

# we can compute the number of non-NA values in our distance metric. Is should be length(chr1$startTAD)
max(ESCchr1$startTAD, ESCchr1$endTAD, na.rm = TRUE)**2 - sum(is.na(ESCchr1metric)) == length(ESCchr1$startTAD)

# compute the eccentricity for each TAD
# we need to include the lower triangular part of ESCchr1metric
# the lower triangular part of chr1metric
for (i in 1:n1) {
        for (j in 1:n1){
                ESCchr1metric[j,i] = ESCchr1metric[i,j]
        }
}

ESCchr1ecc <- matrix(, n1 ,2)

for (i in 1:n1){
        ESCchr1ecc[i,1] = i
        ESCchr1ecc[i,2] = sum(ESCchr1metric[,i], na.rm=TRUE)/sum(!is.na(ESCchr1metric[,i]))
}

ESCchr1ecc <- data.frame(ESCchr1ecc)
names(ESCchr1ecc) = c("TAD", "eccentricity")

# now fill in the rest of the distance metric matrix
# then the diagonal of chr1metric
for (i in 1:n1){ESCchr1metric[i,i] = 0}

# fill in the NA values to be something large? we make them 2 larger than the largest dist value we have
for (i in 1:n1) {
        for (j in 1:n1){
                if (is.na(ESCchr1metric[i,j])){ 
                        ESCchr1metric[i,j] = max(ESCchr1$dist, na.rm = TRUE) + 20
                }
        }
}

# turn chr1metric into a dist metric of class 'dist'
ESCchr1dist = dist(ESCchr1metric)

##########################################################################
# now we do the same as above for Neuronchr1
Neuronchr1col <- which(names(rawNeurondata)=="chr1")
Neuronchr1 <- rawNeurondata[2:length(rawNeurondata[,Neuronchr1col]), (Neuronchr1col+1):(Neuronchr1col+3)]

Neuronchr1[,1] <- as.numeric(Neuronchr1[,1])
Neuronchr1[,2] <- as.numeric(Neuronchr1[,2])
Neuronchr1[,3] <- as.numeric(Neuronchr1[,3])

names(Neuronchr1) = c("startTAD", "endTAD", "interaction")

Neuronchr1 <- Neuronchr1[1:sum(!is.na(ESCchr1$startTAD)), ]

Neuronchr1$dist = 1/Neuronchr1$interaction

# now we want to create a distance metric matrix
n1 = max(Neuronchr1$endTAD, Neuronchr1$startTAD) - (min(Neuronchr1$endTAD, Neuronchr1$startTAD) - 1)
Neuronchr1metric <- matrix(, n1, n1)

for (i in 1:length(Neuronchr1$startTAD)) {
        Neuronchr1metric[Neuronchr1$startTAD[i], Neuronchr1$endTAD[i]] = Neuronchr1$dist[i]
}

# we can compute the number of non-NA values in our distance metric. Is should be length(chr1$startTAD)
max(Neuronchr1$startTAD, Neuronchr1$endTAD, na.rm = TRUE)**2 - sum(is.na(Neuronchr1metric)) == length(Neuronchr1$startTAD)

# compute the eccentricity for each TAD
# we need to include the lower triangular part of ESCchr1metric
# the lower triangular part of chr1metric
for (i in 1:n1) {
        for (j in 1:n1){
                Neuronchr1metric[j,i] = Neuronchr1metric[i,j]
        }
}

Neuronchr1ecc <- matrix(, n1 ,2)

for (i in 1:n1){
        Neuronchr1ecc[i,1] = i
        Neuronchr1ecc[i,2] = sum(Neuronchr1metric[,i], na.rm=TRUE)/sum(!is.na(Neuronchr1metric[,i]))
}

Neuronchr1ecc <- data.frame(Neuronchr1ecc)
names(Neuronchr1ecc) = c("TAD", "eccentricity")

# now fill in the rest of the distance metric matrix
# then the diagonal of chr1metric
for (i in 1:n1){Neuronchr1metric[i,i] = 0}

# fill in the NA values to be something large? we make them 2 larger than the largest dist value we have
for (i in 1:n1) {
        for (j in 1:n1){
                if (is.na(Neuronchr1metric[i,j])){ 
                        Neuronchr1metric[i,j] = max(Neuronchr1$dist, na.rm = TRUE) + 2
                }
        }
}

# turn chr1metric into a dist metric of class 'dist'
Neuronchr1dist = dist(Neuronchr1metric)

####################################################################################################
####################################################################################################
# now use TDA mapper to create a node graph with filter_values = eccentricity
par(mfrow = c(1,2))

ESCchr1.mapper <- mapper1D(distance_matrix = ESCchr1dist,
                           percent_overlap = 75)
ESCchr1.graph <- graph.adjacency(ESCchr1.mapper$adjacency, mode = "undirected")

plot(ESCchr1.graph, layout = layout.auto(ESCchr1.graph))
title(main = "ESCchr1")


Neuronchr1.mapper <- mapper1D(distance_matrix = Neuronchr1dist,
                              filter_values = Neuronchr1ecc,
                              num_intervals = 11)
Neuronchr1.graph <- graph.adjacency(Neuronchr1.mapper$adjacency, mode = "undirected")

plot(Neuronchr1.graph, layout = layout.auto(Neuronchr1.graph))
title(main = "Neuronchr1")
