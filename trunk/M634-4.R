source("common.R")
source("datafile-utils.R")

## Read in data
sys <- read.sys("../../data/Anatomy/ALU/M643-4/CONTRA")
map <- as.matrix(read.map("../../data/Anatomy/ALU/M643-4/CONTRA"))

## Corner analysis
segs <- map.to.segments(map)
## Some curation is required here
segs4 <- segs[[4]][23:1,]
edge.path <- create.path(list(rbind(segs[[1]], segs[[3]], segs4)), close=TRUE)

## Hand pick corners
## c(16, 31)
## c(53, 65)
## c(83, 102)
## Define rim as list of line segments
rim <- list(edge.path[16:31,1:2],
            edge.path[53:65,1:2],
            edge.path[83:94,1:2])

## Find distance along edges
rim.path <- create.path(rim)

cut <- list(edge.path[1:16,1:2],
            edge.path[31:53,1:2],
            edge.path[65:83,1:2],
            edge.path[94:110,1:2])

cut.path <- create.path(cut)

## Experimental code: add a hand-curated tear pair
tears <- list()
tears[[1]]  =
  list(create.path(list(edge.path[42:31,1:2])),
       create.path(list(edge.path[42:53,1:2])))
tears[[2]]  =
  list(create.path(list(edge.path[78:83,1:2])),
       create.path(list(edge.path[78:76,1:2], edge.path[70:67,1:2])))
tears[[3]]  =
  list(create.path(list(edge.path[c(108,1:16),1:2])),
       create.path(list(edge.path[108:94,1:2])))
tears[[4]]  =
  list(create.path(list(edge.path[74:76,1:2])),
       create.path(list(edge.path[74:70,1:2])))
tears[[5]]  =
  list(create.path(list(edge.path[88:89,1:2])),
       create.path(list(edge.path[88:87,1:2])))
tears[[6]]  =
  list(create.path(list(edge.path[66:67,1:2])),
       create.path(list(edge.path[66:65,1:2])))

tearmat <- rbind(c(42,  31, 53),
                 c(78,  67, 83),
                 c(108, 94, 16),
                 c(74,  70, 76),
                 c(88,  87, 89),
                 c(66,  65, 67))
colnames(tearmat) <- c("apex", "end1", "end2")

ifix <- c(16:31, 53:65, 83:94)
