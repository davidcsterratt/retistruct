## Load the good raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "GMB530/R-CONTRA")
r.good <- retistruct.read.dataset(dataset)
## Load the human annotation of tears
r.good <- retistruct.read.markup(r.good)
## Reconstruct
r.good <- retistruct.reconstruct(r.good)

## Load the bad raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "GM182-4/R-CONTRA")
r.bad <- retistruct.read.dataset(dataset)
## Load the human annotation of tears
r.bad <- retistruct.read.markup(r.bad)
## Reconstruct
r.bad <- retistruct.reconstruct(r.bad)

## Plotting
x11(width=6.83, height=6.83/2)
par(mfrow=c(2, 4))

## Good Retina
par(mar=c(2.2, 2.2, 0.5, 0.5),
    mgp=c(1.1, 0.2, 0),
    tcl=-0.15)
lvsLplot(r.good)
panlabel("A")

par(mar=c(0.5, 0.5, 0.5, 0.5))
flatplot(r.good, strain=TRUE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=FALSE, landmarks=FALSE)
panlabel("B")

flatplot(r.good, strain=FALSE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=TRUE)
panlabel("C")

projection(r.good, datapoints=FALSE, datapoint.means=FALSE, datapoint.contours=FALSE)
panlabel("D")

## Bad Retina
par(mar=c(2.2, 2.2, 0.5, 0.5))
lvsLplot(r.bad)
panlabel("E")

par(mar=c(0.5, 0.5, 0.5, 0.5))
flatplot(r.bad, strain=TRUE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=FALSE, landmarks=FALSE)
panlabel("F")

flatplot(r.bad, strain=FALSE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=TRUE)
panlabel("G")

projection(r.bad, datapoints=FALSE, datapoint.means=FALSE, datapoint.contours=FALSE)
panlabel("H")

## Printing
## dev.copy2eps(file="fig2-retistruct-good-bad.eps", width=6.83, height=6.83/2)

