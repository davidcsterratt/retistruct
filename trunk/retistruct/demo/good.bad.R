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
x11(width=12, height=6)
par(mfrow=c(2, 4))

## Good Retina
par(mar=c(2.2, 2.2, 0.5, 0.5),
    mgp=c(1.1, 0.2, 0),
    tcl=-0.15)
plot.l.vs.L(r.good)
mtext("A", adj=-0.2, font=2, line=-0.9)

par(mar=c(0.5, 0.5, 0.5, 0.5))
plot.flat(r.good, strain=TRUE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=FALSE, landmarks=FALSE)
mtext("B", adj=0, font=2, line=-0.9)
plot.flat(r.good, strain=FALSE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=TRUE)
mtext("C", adj=0, font=2, line=-0.9)
plot.polar(r.good, datapoints=FALSE, datapoint.means=FALSE, datapoint.contours=FALSE)
mtext("D", adj=0, font=2, line=-0.9)

## Bad Retina
par(mar=c(2.2, 2.2, 0.5, 0.5))
plot.l.vs.L(r.bad)
mtext("E", adj=-0.2, font=2, line=-0.9)

par(mar=c(0.5, 0.5, 0.5, 0.5))
plot.flat(r.bad, strain=TRUE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=FALSE, landmarks=FALSE)
mtext("F", adj=0, font=2, line=-0.9)
plot.flat(r.bad, strain=FALSE, mesh=FALSE, stitch=FALSE, markup=FALSE, datapoints=FALSE, grid=TRUE)
mtext("G", adj=0, font=2, line=-0.9)
plot.polar(r.bad, datapoints=FALSE, datapoint.means=FALSE, datapoint.contours=FALSE)
mtext("H", adj=0, font=2, line=-0.9)

## Printing
## dev.copy2pdf(file="retistruct-good-bad.pdf", width=6.83, height=6.83/2)
