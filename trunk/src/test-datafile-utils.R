source("datafile-utils.R")

datasets <- c("../../../data/Anatomy/marked-up-retinae-2010-03-24/gm119-5-adult-C57BL6",
              "../../../data/Anatomy/marked-up-retinae-2010-03-24/gm143-2-P2-C57BL6/",
              "../../../data/Anatomy/marked-up-retinae-2010-03-24/gm257-1-P8-C57BL6/",
              "../../../data/Anatomy/ALU/M642-3/CONTRA")

##../ 	./ 	CONFIG/ 	M641-5/ 	M641-6/ 	M642-1/ 	M642-2/
##M642-3/ 	M642-4/ 	M642-5/ 	M642-6/ 	M642-7/ 	M643-1/ 	M643-2/
##M643-3/ 	M643-4/ 	M643-5/ 	M643-6/ 	m643.pdf

for (dataset in datasets) {
  map <- read.map(dataset)
  ## sys <- read.sys(dataset)
  plot.map(map, TRUE)

  segs <- map.to.segments(map)
  P <- segments.to.outline(segs)
  lines(P[,1], P[,2])
  text(P[,1], P[,2]-100, 1:nrow(P))
  readline("Press <Enter> to continue")
}
