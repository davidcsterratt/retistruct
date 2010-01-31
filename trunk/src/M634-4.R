source("fold-sphere2.R")
source("common.R")
source("datafile-utils.R")

## Read in data
sys <- read.sys("data/Anatomy/ALU/M643-4/CONTRA")
map <- as.matrix(read.map("data/Anatomy/ALU/M643-4/CONTRA"))

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

## Actual munging routine
P <- edge.path[-nrow(edge.path),1:2]

s <- stitch.retina(P, tearmat)
plot.stitch(s)

t <- make.triangulation(s)
with(t, trimesh(T, P, col="black"))

m <- merge.points(t, s)
## Plot stiched retina in 2D (messy)
## trimesh(Tt, Pt, col="black")

## Plotting
plot(P)
with(s, plot.outline(P, gb))

p <- project.to.sphere(m, t, phi0=50*pi/180)

## Initial plot in 3D space
plot.retina(p$phi, p$lambda, p$R, m$Tt, m$Rsett)

## In order to plot the points, we need to know in which triangle they
## lie and where in that triangle they are. The first job is to create
## a triangulation, and then we can use tsearch to find the identity of the
## triangles and the location in the triangles.
P.red   <- cbind(na.omit(sys[,"XRED"]), na.omit(sys[, "YRED"]))
P.green <- cbind(na.omit(sys[,"XGREEN"]), na.omit(sys[, "YGREEN"]))
cb.red   <- with(t, tsearchn(P, T, P.red))
cb.green <- with(t, tsearchn(P, T, P.green))

## Attempt to solve (or refine?) mapping
# r <- solve.mapping(p, m, t, s, E0.A=0, dt=2E-6, nstep=5000, Rexp=1, verbose=FALSE)
# 3304 after 1000 steps
# 3107 after 4000 steps; there is some jumping around if the tolerance isn't small enough
# 2831 after 4000 steps with lower tolerance 0.002

#r = solve.mapping.momentum(p, m, t, s, E0.A=0, dt=1E-4, nstep=1000, Rexp=1, verbose=FALSE)
# 3439 after 1000 steps
# 3039 after 2000 steps; optimise mapping doesn't work after this
# 2884 after 4000 steps
# 2991 after 4000 steps with lower tolerance



## p1 <- p
## p1$phi <- r$phi
## p1$lambda <- r$lambda

## r1 <- optimise.mapping.nlm(p1, m, t, s, E0.A=0, iterlim=10000)
