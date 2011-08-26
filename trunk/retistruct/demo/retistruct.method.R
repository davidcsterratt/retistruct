## Set up a 3x3 grid for plotting
par(mfrow=c(3, 3))

## Load the raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "GM114-4-RC")
o <- retistruct.read.dataset(dataset)

## Load the human annotation of tears
o <- retistruct.read.markup(o)

## Plot of raw data. Axes are reversed to improve comparison with
## polar plot later
xlim <- range(o$P[,1])[c(2,1)]
ylim <- range(o$P[,2])[c(2,1)]
plot.flat(o, markup=FALSE,
          xlim=xlim, ylim=ylim)
mtext("A", adj=0, font=2, line=-0.9)

## Plot the annotation
plot.flat(o, datapoints=FALSE, landmarks=FALSE,
          xlim=xlim, ylim=ylim)
mtext("B", adj=0, font=2, line=-0.9)

## Set up fixed point
o$lambda0 <- 0

## Initial triangulation (with 500 points)
n <- 500
t <- triangulate.outline(o, n=n)

## Stitching
s <- stitch.outline(t)

## Plot triangulation and stitching
plot.flat(s, datapoints=FALSE, landmarks=FALSE, markup=FALSE,
          xlim=xlim, ylim=ylim)
mtext("C", adj=0, font=2, line=-0.9)

## Triangulate again, to take into account points added by stitching
m <- triangulate.outline(s, n=n,
                         suppress.external.steiner=TRUE)

## Merge the points that have been stitched
m <- merge.points.edges(m)

## Make a rough projection to a sphere
m <- project.to.sphere(m)

## Plot the initial gridlines in 2D
par(mfg=c(3, 3))
plot.flat(m, grid=TRUE, 
          datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE,
          stitch=FALSE, strain=TRUE,
          xlim=xlim, ylim=ylim)
mtext("Dii", adj=0, font=2, line=-0.9)

## Plot the intial projection in 3D
par(mfg=c(2, 3))
plot.new()
mtext("Di", adj=0, font=2, line=-0.9)

plot.sphere.spherical(m$phi, m$lambda, m$R, m$Tt, m$Rsett)
plot.outline.spherical(m$phi, m$lambda, m$R, m$gb, m$ht)
rgl.viewpoint(zoom=0.7)
rgl.postscript("initial-projection.pdf", "pdf")

## Optimise mapping - this takes a few minutes
alpha <- 8
x0 <- 0.6
r <- solve.mapping.cart(m, alpha=alpha, x0=x0, nu=1, dtmax=500, maxmove=1E3,
                        tol=1e-5, plot.3d=FALSE)
r <- optimise.mapping(r, alpha=alpha, x0=x0, nu=0.5, 
                      plot.3d=FALSE)

## Plot the final projection in 3D and on the grid
par(mfg=c(2, 2))
plot.new()
mtext("Ei", adj=0, font=2, line=-0.9)

rgl.postscript("final-projection.pdf", "pdf")

par(mfg=c(3, 2))
plot.flat(r, grid=TRUE, 
          datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE,
          stitch=FALSE, strain=TRUE,
          xlim=xlim, ylim=ylim)
mtext("Eii", adj=0, font=2, line=-0.9)

## Infer locations of datapoints in spherical coordinates
r <- ReconstructedDataset(r)
r <- RetinalReconstructedDataset(r)

## Plot data in polar coordinates and flattend retina
par(mfg=c(2, 1))
plot.polar(r, datapoints=TRUE, landmarks=TRUE)
mtext("Fi", adj=0, font=2, line=-0.9)

par(mfg=c(3, 1))
plot.flat(r, grid=TRUE, 
          datapoints=TRUE, landmarks=TRUE, mesh=FALSE, markup=FALSE,
          stitch=FALSE,
          xlim=xlim, ylim=ylim)
mtext("Fii", adj=0, font=2, line=-0.9)

## Save to PDF
dev.copy2pdf(file="retistruct-method.pdf", width=6.83, height=6.83)
