## Load the raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "smi32")
o <- retistruct.read.dataset(dataset)

## Load the human annotation of tears
o <- retistruct.read.markup(o)

## Initial plot
par(mar=c(0.1,0.1,0.1,0.1))
plot.flat(o)

## Reconstruct
r <- retistruct.reconstruct(o)

## Plot with gridlines
options(outline.col="yellow")
options(grid.maj.col="yellow")
options(grid.min.col="white")
plot.flat(r, mesh=FALSE, stitch=FALSE, markup=FALSE)

## To print for PloS Biol.
## dev.copy2eps(file="KeelEtal11.eps", width=6.83, height=6.83)
