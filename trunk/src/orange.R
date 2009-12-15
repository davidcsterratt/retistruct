library("pixmap")

im <- read.pnm("../data/20091214191056-small.ppm")
# plot(im)

## Read in data
P.raw <- read.table("../data/20091214191056-small.txt")
P <- cbind(P.raw[,1], 520-P.raw[,2])

tearmat <- rbind(c(5,   1, 11),
                 c(21, 16, 32),
                 c(29, 25, 31),
                 c(52, 39, 65))
colnames(tearmat) <- c("apex", "end1", "end2")


