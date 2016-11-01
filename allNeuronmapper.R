# this script creates a node graph for chr1, chr2, chr3, chr4 of the Neuron data
# all setting in the mapper1D function are set to default values
# filter_values = a length n vector of real numbers (default?)
# num_intervals = 10
# percent_overlap = 50
# num_bins_when_clustering = 10

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


rawNeurondata <- read.csv("~/Dropbox/DS/wTrush/rawNeurondata.csv", header=TRUE, stringsAsFactors = FALSE)
dim(rawNeurondata)
str(rawNeurondata)
# we can see this is a 415 x 16 dataframe
# start by isolating the relevant data for chr1
chr1col <- which(names(rawNeurondata)=="chr1")
chr1 <- rawNeurondata[2:length(rawESCdata[,chr1col]), (chr1col+1):(chr1col+3)]

chr1[,1] <- as.numeric(chr1[,1])
chr1[,2] <- as.numeric(chr1[,2])
chr1[,3] <- as.numeric(chr1[,3])

names(chr1) <- c("startTAD", "endTAD", "interaction")

# this isn't a problem for chr1 because it has the most TAD information. However, we want to make sure we ignore any NA
# values that might be at the end of these columns. 
chr1 <- chr1[1:sum(!is.na(chr1$startTAD)), ]

# next we create a new column which defines 'distance' between given TADs
# here we use 1/interaction
chr1$dist = 1/chr1$interaction

# note the max and min values
max(chr1$dist)
min(chr1$dist)

# now we want to create a distance metric
max(chr1$endTAD, chr1$startTAD)
min(chr1$endTAD, chr1$startTAD)
# given the max/min values of the TADs available, we will want a 415x415 matrix

n1 = max(chr1$endTAD, chr1$startTAD) - (min(chr1$endTAD, chr1$startTAD) - 1)
chr1metric <- matrix(, n1, n1)

for (i in 1:length(chr1$startTAD)) {
        chr1metric[chr1$startTAD[i], chr1$endTAD[i]] = chr1$dist[i]
}

# check out chr1metric
head(chr1metric)

# we can compute the number of non-NA values in our distance metric. Is should be length(chr1$startTAD)
max(chr1$startTAD, chr1$endTAD, na.rm = TRUE)**2 - sum(is.na(chr1metric)) == length(chr1$startTAD)

# also note, this is an upper triangular matrix. we can fill in the diagonal and lower triangular part
# to get the `full' distance metric on chr1

# 1st the lower triangular part
for (i in 1:n1) {
        for (j in 1:n1){
                chr1metric[j,i] = chr1metric[i,j]
        }
}
# then the diagonal
for (i in 1:n1){ chr1metric[i,i] = 0}

# now fill in the NA values to be something large? we make them 2 larger than the largest dist value we have
for (i in 1:n1) {
        for (j in 1:n1){
                if (is.na(chr1metric[i,j])){ 
                        chr1metric[i,j] = max(chr1$dist, na.rm = TRUE) + 200
                }
        }
}

# now create the distance metric
# dist creates an object in the class `dist' in R
chr1dist = dist(chr1metric)

Neuronchr1.mapper <- mapper1D(distance_matrix = chr1dist)
Neuronchr1.graph <- graph.adjacency(Neuronchr1.mapper$adjacency, mode = "undirected")

par(mfrow = c(2,2))
plot(Neuronchr1.graph, layout = layout.auto(Neuronchr1.graph))
title(main = "Neuron chr1")

########################################################        
# Ok. we get a node graph with the default values. 
# Let's see how different the node graphs look for all the neurons
########################################################        

########################################################        
##############         chr2             ################
########################################################        
# start by isolating the relevant data for chr2. note the code below is more adaptable/flexible
chr2col <- which(names(rawNeurondata)=="chr2")
chr2 <- rawNeurondata[2:length(rawNeurondata[,chr2col]), (chr2col+1):(chr2col+3)]

chr2[,1] <- as.numeric(chr2[,1])
chr2[,2] <- as.numeric(chr2[,2])
chr2[,3] <- as.numeric(chr2[,3])

names(chr2) <- c("startTAD", "endTAD", "interaction")

# because we have less data for chr2, there are lots of NA values at the end of the columns for chr2.
tail(chr2)
# we want to remove these NA values. So let's find out how many values we have 
sum(!is.na(chr2$startTAD))
# so there are only 359 TAD interaction measurements. let's isolate those
chr2 <- chr2[1:sum(!is.na(chr2$startTAD)), ]

chr2$dist = 1/chr2$interaction

# note the max and min values
max(chr2$dist)
min(chr2$dist)

# now we want to create a distance metric
max(chr2$endTAD, chr2$startTAD)
min(chr2$endTAD, chr2$startTAD)
# given the above results of the max/min values of the TADs available, we will want a 373x373 matrix

n2 = max(chr2$endTAD, chr2$startTAD) - min(chr2$endTAD, chr2$startTAD) + 1
chr2metric <- matrix(, n2, n2)

for (i in 1:n2) {
        chr2metric[chr2$startTAD[i], chr2$endTAD[i]] = chr2$dist[i]
}

# check out chr2metric
head(chr2metric)
tail(chr2metric)

# check how many values we do have. the computed value below should be length(chr2$startTAD)
max(chr2$startTAD, chr2$endTAD)**2 - sum(is.na(chr2metric)) == length(chr2$startTAD)

# now fill in the diagonal and lower triangular part
# 1st the lower triangular part
for (i in 1:n2) {
        for (j in 1:n2){
                chr2metric[j,i] = chr2metric[i,j]
        }
}
# then the diagonal
for (i in 1:n2){ chr2metric[i,i] = 0}

# now fill in the NA values to be something large? 
for (i in 1:n2) {
        for (j in 1:n2){
                if (is.na(chr2metric[i,j])){ 
                        chr2metric[i,j] = max(chr2$dist) + 2
                }
        }
}

# now create the distance metric
chr2dist = dist(chr2metric)

Neuronchr2.mapper <- mapper1D(distance_matrix = chr2dist)
Neuronchr2.graph <- graph.adjacency(Neuronchr2.mapper$adjacency, mode = "undirected")
plot(Neuronchr2.graph, layout = layout.auto(Neuronchr2.graph))
title(main = "Neuron chr2")

########################################################        
##############         chr3             ################
########################################################        
chr3col <- which(names(rawNeurondata)=="chr3")
chr3 <- rawNeurondata[2:length(rawNeurondata[,chr3col]), (chr3col+1):(chr3col+3)]

chr3[,1] <- as.numeric(chr3[,1])
chr3[,2] <- as.numeric(chr3[,2])
chr3[,3] <- as.numeric(chr3[,3])

names(chr3) <- c("startTAD", "endTAD", "interaction")

# we also have less data for chr3 and we want to remove these NA values. 
chr3 <- chr3[1:sum(!is.na(chr3$startTAD)), ]

chr3$dist = 1/chr3$interaction

# note the max and min values
max(chr3$dist)
min(chr3$dist)

# now we want to create a distance metric
max(chr3$startTAD, chr3$endTAD)
min(chr3$startTAD, chr3$endTAD)

n3 = max(chr3$endTAD, chr3$startTAD) - (min(chr3$endTAD, chr3$startTAD) - 1)
chr3metric <- matrix(, n3, n3)

for (i in 1:n3) {
        chr3metric[chr3$startTAD[i], chr3$endTAD[i]] = chr3$dist[i]
}

# check out chr3dist. 
head(chr3metric)
tail(chr3metric)

# double-check how many values we have so far. the computed value below should be length(chr3$startTAD)
max(chr3$startTAD, chr3$endTAD)**2 - sum(is.na(chr3metric)) == length(chr3$startTAD)

# now fill in the diagonal and lower triangular part
# 1st the lower triangular part
for (i in 1:n3) {
        for (j in 1:n3){
                chr3metric[j,i] = chr3metric[i,j]
        }
}
# then the diagonal
for (i in 1:n3){ chr3metric[i,i] = 0}

# now fill in the NA values to be something large? 
for (i in 1:n3) {
        for (j in 1:n3){
                if (is.na(chr3metric[i,j])){ 
                        chr3metric[i,j] = max(chr3$dist) + 2
                }
        }
}

# distance metric for chr3
chr3dist = dist(chr3metric)

Neuronchr3.mapper <- mapper1D(distance_matrix = chr3dist)
Neuronchr3.graph <- graph.adjacency(Neuronchr3.mapper$adjacency, mode = "undirected")
plot(Neuronchr3.graph, layout = layout.auto(Neuronchr3.graph))
title(main = "Neuron chr3")

########################################################        
##############         chr4             ################
########################################################        
chr4col <- which(names(rawNeurondata)=="chr4")
chr4 <- rawNeurondata[2:length(rawNeurondata[,chr4col]), (chr4col+1):(chr4col+3)]

chr4[,1] <- as.numeric(chr4[,1])
chr4[,2] <- as.numeric(chr4[,2])
chr4[,3] <- as.numeric(chr4[,3])

names(chr4) = c("startTAD", "endTAD", "interaction")

# we also have less data for chr4 and we want to remove these NA values. 
chr4 <- chr4[1:sum(!is.na(chr4$startTAD)), ]

chr4$dist = 1/chr4$interaction

# note the max and min values
max(chr4$dist)
min(chr4$dist)

# now we want to create a distance metric
max(chr4$startTAD, chr4$endTAD)
min(chr4$startTAD, chr4$endTAD)

n4 = max(chr4$endTAD, chr4$startTAD) - (min(chr4$endTAD, chr4$startTAD) - 1)
chr4metric <- matrix(, n4, n4)

for (i in 1:n4) {
        chr4metric[chr4$startTAD[i], chr4$endTAD[i]] = chr4$dist[i]
}

# check out chr4dist. 
head(chr4metric)
tail(chr4metric)

# double-check how many values we have so far. the computed value below should be length(chr4$startTAD)
max(chr4$startTAD, chr4$endTAD)**2 - sum(is.na(chr4metric)) == length(chr4$startTAD)

# now fill in the diagonal and lower triangular part
# 1st the lower triangular part
for (i in 1:n4) {
        for (j in 1:n4){
                chr4metric[j,i] = chr4metric[i,j]
        }
}
# then the diagonal
for (i in 1:n4){ chr4metric[i,i] = 0}

# now fill in the NA values to be something large? 
for (i in 1:n4) {
        for (j in 1:n4){
                if (is.na(chr4metric[i,j])){ 
                        chr4metric[i,j] = max(chr4$dist) + 2
                }
        }
}

# distance metric for chr4
chr4dist = dist(chr4metric)

Neuronchr4.mapper <- mapper1D(distance_matrix = chr4dist)
Neuronchr4.graph <- graph.adjacency(Neuronchr4.mapper$adjacency, mode = "undirected")
plot(Neuronchr4.graph, layout = layout.auto(Neuronchr4.graph))
title(main = "Neuron chr4")

########################################################        
##############    Plotted individually      ############
########################################################        

par(mfrow=c(1,1))
plot(Neuronchr1.graph, layout = layout.auto(Neuronchr1.graph))
title(main = "Neuron chr1") 

par(mfrow=c(1,1))
plot(Neuronchr2.graph, layout = layout.auto(Neuronchr2.graph))
title(main = "Neuron chr2") 

par(mfrow=c(1,1))
plot(Neuronchr3.graph, layout = layout.auto(Neuronchr3.graph))
title(main = "Neuron chr3") 

par(mfrow=c(1,1))
plot(Neuronchr4.graph, layout = layout.auto(Neuronchr4.graph))
title(main = "Neuron chr4") 


########################################################        
#########       ESC/Neuron Plotted together      #######
########################################################        
# note for this segment of code we need all ESC and Neuron mapper results

par(mfrow=c(1,2))
plot(ESCchr1.graph, layout = layout.auto(ESCchr1.graph))
title(main = "ESC chr1")
plot(Neuronchr1.graph, layout = layout.auto(Neuronchr1.graph))
title(main = "Neuron chr1")

par(mfrow=c(1,2))
plot(ESCchr2.graph, layout = layout.auto(ESCchr2.graph))
title(main = "ESC chr2")
plot(Neuronchr2.graph, layout = layout.auto(Neuronchr2.graph))
title(main = "Neuron chr2")

par(mfrow=c(1,2))
plot(ESCchr3.graph, layout = layout.auto(ESCchr3.graph))
title(main = "ESC chr3")
plot(Neuronchr3.graph, layout = layout.auto(Neuronchr3.graph))
title(main = "Neuron chr3")

par(mfrow=c(1,2))
plot(ESCchr4.graph, layout = layout.auto(ESCchr4.graph))
title(main = "ESC chr4")
plot(Neuronchr4.graph, layout = layout.auto(Neuronchr4.graph))
title(main = "Neuron chr4")

