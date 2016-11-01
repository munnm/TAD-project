library(ggplot2)

rawESCdata <- read.csv("~/Dropbox/DS/wTrush/rawESCdata.csv", header=TRUE, stringsAsFactors = FALSE)

chr1 <- rawESCdata[2:length(rawESCdata[,2]), 2:4]

chr1$X <- as.numeric(chr1$X)
chr1$X.1 <- as.numeric(chr1$X.1)
chr1$X.2 <- as.numeric(chr1$X.2)

names(chr1) = c("startTAD", "endTAD", "interaction")

chr1 <- chr1[1:sum(!is.na(chr1$startTAD)), ]

chr1$dist = 1/chr1$interaction

n1 = max(chr1$endTAD, chr1$startTAD) - (min(chr1$endTAD, chr1$startTAD) - 1)
chr1metric <- matrix(, n1, n1)

for (i in 1:length(chr1$startTAD)) {
        chr1metric[chr1$startTAD[i], chr1$endTAD[i]] = chr1$dist[i]
}

# we can compute the number of non-NA values in our distance metric. Is should be length(chr1$startTAD)
max(chr1$startTAD, chr1$endTAD, na.rm = TRUE)**2 - sum(is.na(chr1metric)) == length(chr1$startTAD)

# 1st the lower triangular part
for (i in 1:n1) {
        for (j in 1:n1){
                chr1metric[j,i] = chr1metric[i,j]
        }
}
# we don't need to do the diagonal since each diagonal value is 0 and affects the number of entries

chr1ecc <- matrix(, 394,2)

for (i in 1:394){
        chr1ecc[i,1] = i
        chr1ecc[i,2] = sum(chr1metric[,i], na.rm=TRUE)/sum(!is.na(chr1metric[,i]))
}

chr1ecc <- data.frame(chr1ecc)
names(chr1ecc) = c("TAD", "eccentricity")

head(chr1ecc)

# see if any TADs lack a value
sum(is.na(chr1ecc$eccentricity))

par(mfrow = c(1,2))
g <- ggplot(chr1ecc, aes(TAD, log(eccentricity)))
p <- g + geom_point(alpha = 1/2, size = 2) + geom_smooth(method = "lm")  + ggtitle("chr1 eccentricity")
print(p)
########################################################        
##############         chr2             ################
########################################################        
# start by isolating the relevant data for chr2. note the code below is more adaptable/flexible
chr2col <- which(names(rawESCdata)=="chr2")
chr2 <- rawESCdata[2:length(rawESCdata[,chr2col]), (chr2col+1):(chr2col+3)]

chr2[,1] <- as.numeric(chr2[,1])
chr2[,2] <- as.numeric(chr2[,2])
chr2[,3] <- as.numeric(chr2[,3])

names(chr2) <- c("startTAD", "endTAD", "interaction")

# because we have less data for chr2, there are lots of NA values at the end of the columns for chr2.
tail(chr2)
# we want to remove these NA values. So let's find out how many values we have 
sum(!is.na(chr2$startTAD))
# so there are only 372 TAD interaction measurements. let's isolate those
chr2 <- chr2[1:sum(!is.na(chr2$startTAD)), ]

chr2$dist = 1/chr2$interaction

n2 = max(chr2$endTAD, chr2$startTAD) - min(chr2$endTAD, chr2$startTAD) + 1

chr2metric <- matrix(, n2, n2)

for (i in 1:n2) {
        chr2metric[chr2$startTAD[i], chr2$endTAD[i]] = chr2$dist[i]
}

# check how many values we do have. the computed value below should be length(chr2$startTAD)
max(chr2$startTAD, chr2$endTAD)**2 - sum(is.na(chr2metric)) == length(chr2$startTAD)

# now fill in the diagonal and lower triangular part
# 1st the lower triangular part
for (i in 1:n2) {
        for (j in 1:n2){
                chr2metric[j,i] = chr2metric[i,j]
        }
}

chr2ecc <- matrix(, n2,2)

for (i in 1:n2){
        chr2ecc[i,1] = i
        chr2ecc[i,2] = sum(chr2metric[,i], na.rm=TRUE)/sum(!is.na(chr2metric[,i]))
}


chr2ecc <- data.frame(chr2ecc)
names(chr2ecc) = c("TAD", "eccentricity")

head(chr2ecc)

# see if any TADs lack a value
sum(is.na(chr2ecc$eccentricity))


par(mfrow = c(1,2))
g <- ggplot(chr2ecc, aes(TAD, log(eccentricity)))
p <- g + geom_point(alpha = 1/2, size = 2) + geom_smooth(method = "lm") + ggtitle("chr2 eccentricity")
print(p)

########################################################        
##############         chr3             ################
########################################################        
# start by isolating the relevant data for chr2. note the code below is more adaptable/flexible
chr3col <- which(names(rawESCdata)=="chr3")
chr3 <- rawESCdata[2:length(rawESCdata[,chr3col]), (chr3col+1):(chr3col+3)]

chr3[,1] <- as.numeric(chr3[,1])
chr3[,2] <- as.numeric(chr3[,2])
chr3[,3] <- as.numeric(chr3[,3])

names(chr3) <- c("startTAD", "endTAD", "interaction")

chr3 <- chr3[1:sum(!is.na(chr3$startTAD)), ]

chr3$dist = 1/chr3$interaction

n3 = max(chr3$endTAD, chr3$startTAD) - min(chr3$endTAD, chr3$startTAD) + 1

chr3metric <- matrix(, n3, n3)

for (i in 1:n3) {
        chr3metric[chr3$startTAD[i], chr3$endTAD[i]] = chr3$dist[i]
}

# check how many values we do have. the computed value below should be length(chr2$startTAD)
max(chr3$startTAD, chr3$endTAD)**2 - sum(is.na(chr3metric)) == length(chr3$startTAD)

# now fill in the diagonal and lower triangular part
# 1st the lower triangular part
for (i in 1:n3) {
        for (j in 1:n3){
                chr3metric[j,i] = chr3metric[i,j]
        }
}

chr3ecc <- matrix(, n3,2)

for (i in 1:n3){
        chr3ecc[i,1] = i
        chr3ecc[i,2] = sum(chr3metric[,i], na.rm=TRUE)/sum(!is.na(chr3metric[,i]))
}


chr3ecc <- data.frame(chr3ecc)
names(chr3ecc) = c("TAD", "eccentricity")

head(chr3ecc)

# see if any TADs lack a value
sum(is.na(chr3ecc$eccentricity))

par(mfrow = c(1,2))
g <- ggplot(chr3ecc, aes(TAD, log(eccentricity)))
p <- g + geom_point(alpha = 1/2, size = 2) + geom_smooth(method = "lm") + ggtitle("chr3 eccentricity")
print(p)

########################################################        
##############         chr4             ################
########################################################        
# start by isolating the relevant data for chr2. note the code below is more adaptable/flexible
chr4col <- which(names(rawESCdata)=="chr4")
chr4 <- rawESCdata[2:length(rawESCdata[,chr4col]), (chr4col+1):(chr4col+3)]

chr4[,1] <- as.numeric(chr4[,1])
chr4[,2] <- as.numeric(chr4[,2])
chr4[,3] <- as.numeric(chr4[,3])

names(chr4) <- c("startTAD", "endTAD", "interaction")

chr4 <- chr4[1:sum(!is.na(chr4$startTAD)), ]

chr4$dist = 1/chr4$interaction

n4 = max(chr4$endTAD, chr4$startTAD) - min(chr4$endTAD, chr4$startTAD) + 1

chr4metric <- matrix(, n4, n4)

for (i in 1:n4) {
        chr4metric[chr4$startTAD[i], chr4$endTAD[i]] = chr4$dist[i]
}

# check how many values we do have. the computed value below should be length(chr2$startTAD)
max(chr4$startTAD, chr4$endTAD)**2 - sum(is.na(chr4metric)) == length(chr4$startTAD)

# now fill in the diagonal and lower triangular part
# 1st the lower triangular part
for (i in 1:n4) {
        for (j in 1:n4){
                chr4metric[j,i] = chr4metric[i,j]
        }
}

chr4ecc <- matrix(, n4,2)

for (i in 1:n4){
        chr4ecc[i,1] = i
        chr4ecc[i,2] = sum(chr4metric[,i], na.rm=TRUE)/sum(!is.na(chr4metric[,i]))
}


chr4ecc <- data.frame(chr4ecc)
names(chr4ecc) = c("TAD", "eccentricity")

head(chr4ecc)

# see if any TADs lack a value
sum(is.na(chr4ecc$eccentricity))

par(mfrow = c(1,2))
g <- ggplot(chr4ecc, aes(TAD, log(eccentricity)))
p <- g + geom_point(alpha = 1/2, size = 2) + geom_smooth(method = "lm") + ggtitle("chr 4 eccentricity")
print(p)

