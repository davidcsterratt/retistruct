## Read in data
sys <- read.sys("../../data/Anatomy/ALU/M643-4/CONTRA")
map <- as.matrix(read.map("../../data/Anatomy/ALU/M643-4/CONTRA"))

## Corner analysis
segs <- map.to.segments(map)
## Some curation is required here
segs4 <- segs[[4]][23:1,]
edge.path <- create.path(list(segs[[1]], segs[[3]], segs4), close=TRUE)
