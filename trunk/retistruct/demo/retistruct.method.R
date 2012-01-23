## Set up a 3x3 grid for plotting
par(mfrow=c(3, 3))
par(mar=c(0.5, 0.5, 0.5, 0.5))

## Load the raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "GM114-4-RC")
o <- retistruct.read.dataset(dataset)

## Load the human annotation of tears
o <- retistruct.read.markup(o)

## Make this a left eye to help with orientatio of points
o$side="Left"

## Plot of raw data. Axes are reversed to improve comparison with
## polar plot later
plot.flat(o, markup=FALSE)
mtext("A", adj=0, font=2, line=-0.9)

## Plot the annotation
plot.flat(o, datapoints=FALSE, landmarks=FALSE)
mtext("B", adj=0, font=2, line=-0.9)

## Set up fixed point
o$lambda0 <- 0

## Initial triangulation (with 500 points)
n <- 500
t <- TriangulatedOutline(o, n=n)

## Stitching
s <- StitchedOutline(t)

## Plot triangulation and stitching
plot.flat(s, datapoints=FALSE, landmarks=FALSE, markup=FALSE)
mtext("C", adj=0, font=2, line=-0.9)

## Triangulate again, to take into account points added by stitching
m <- TriangulatedOutline(s, n=n,
                         suppress.external.steiner=TRUE)

## Merge the points that have been stitched
m <- merge.points.edges(m)

## Make a rough projection to a sphere
m <- project.to.sphere(m)

## Plot the initial gridlines in 2D
par(mfg=c(3, 3))
plot.flat(m, grid=TRUE, 
          datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE,
          stitch=FALSE, strain=TRUE)
mtext("Dii", adj=0, font=2, line=-0.9)

## Plot the intial projection in 3D
par(mfg=c(2, 3))
plot.new()
mtext("Di", adj=0, font=2, line=-0.9)

plot.spherical(m, strain=TRUE, datapoints=FALSE)
rgl.viewpoint(zoom=0.7)
rgl.postscript("initial-projection.pdf", "pdf")

## Optimise mapping - this takes a few minutes
alpha <- 8
x0 <- 0.5
r <- solve.mapping.cart(m, alpha=0, x0=0, nu=1,
                        dtmax=500, maxmove=1E2, tol=2e-7,
                          plot.3d=FALSE)
r <- solve.mapping.cart(r, alpha=alpha, x0=x0, nu=1,
                        dtmax=500, maxmove=1E2, tol=1e-6,
                        plot.3d=FALSE)
r <- optimise.mapping(r, alpha=alpha, x0=x0, nu=0.5, 
                      plot.3d=FALSE)

## Plot the final projection in 3D and on the grid
par(mfg=c(2, 2))
plot.new()
mtext("Ei", adj=0, font=2, line=-0.9)

plot.spherical(r, strain=TRUE, datapoints=FALSE)
rgl.postscript("final-projection.pdf", "pdf")

par(mfg=c(3, 2))
plot.flat(r, grid=TRUE, 
          datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE,
          stitch=FALSE, strain=TRUE)
mtext("Eii", adj=0, font=2, line=-0.9)

## Infer locations of datapoints in spherical coordinates
r <- ReconstructedDataset(r)
r <- RetinalReconstructedDataset(r)

## Plot data in polar coordinates and flattend retina
par(mfg=c(2, 1))
plot.polar(r, datapoints=TRUE, landmarks=TRUE, datapoint.contours=FALSE)
mtext("Fi", adj=0, font=2, line=-0.9)

par(mfg=c(3, 1))
plot.flat(r, grid=TRUE, 
          datapoints=TRUE, landmarks=TRUE, mesh=FALSE, markup=FALSE,
          stitch=FALSE)
mtext("Fii", adj=0, font=2, line=-0.9)

## Save to PDF
dev.copy2pdf(file="retistruct-method.pdf", width=6.83, height=6.83)
