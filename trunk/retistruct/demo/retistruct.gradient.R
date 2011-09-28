## Load the raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "ChenEtal95")
o <- retistruct.read.dataset(dataset)

## Load the human annotation of tears
o <- retistruct.read.markup(o)

## Initial plot
## plot.flat(o)

## Reconstruct
r <- retistruct.reconstruct(o)

## Plot with gridlines
par(mar=c(0.1,0.1,0.1,0.1))
par(mfrow=c(1, 2))
plot.flat(r, mesh=FALSE, stitch=FALSE, markup=FALSE)
plot.polar(r, mesh=FALSE)

## To print for PloS Biol.
## dev.copy2eps(file="ChenEtal95.eps", width=6.83, height=6.83/2)
