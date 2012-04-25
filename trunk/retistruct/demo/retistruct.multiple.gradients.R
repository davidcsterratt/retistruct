par(mar=c(0.2,0.2,0.2,0.2))
par(mfrow=c(2, 4))

## Load the raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "MarcEtal96EphA")
o <- retistruct.read.dataset(dataset)
## Load the human annotation of tears
o <- retistruct.read.markup(o)
## Reconstruct
r.EphA <- retistruct.reconstruct(o)

## Plot with gridlines
plot.flat(r.EphA, mesh=FALSE, stitch=FALSE, markup=FALSE)
mtext("A", adj=0, font=2, line=-0.9)
polarplot(r.EphA, mesh=FALSE)
mtext("B", adj=0, font=2, line=-0.9)

## Load the raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "MarcEtal96ephrinA")
o <- retistruct.read.dataset(dataset)
## Load the human annotation of tears
o <- retistruct.read.markup(o)
## Reconstruct
r.ephrinA <- retistruct.reconstruct(o)

## Plot with gridlines
plot.flat(r.ephrinA, mesh=FALSE, stitch=FALSE, markup=FALSE)
mtext("C", adj=0, font=2, line=-0.9)
polarplot(r.ephrinA, mesh=FALSE)
mtext("D", adj=0, font=2, line=-0.9)

## Load the raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "MarcEtal96EphB")
o <- retistruct.read.dataset(dataset)
## Load the human annotation of tears
o <- retistruct.read.markup(o)
## Reconstruct
r.EphB <- retistruct.reconstruct(o)

## Plot with gridlines
plot.flat(r.EphB, mesh=FALSE, stitch=FALSE, markup=FALSE)
mtext("E", adj=0, font=2, line=-0.9)
polarplot(r.EphB, mesh=FALSE)
mtext("F", adj=0, font=2, line=-0.9)

## Load the raw data
dataset <- file.path(system.file(package = "retistruct"), "extdata", "MarcEtal96ephrinB")
o <- retistruct.read.dataset(dataset)
## Load the human annotation of tears
o <- retistruct.read.markup(o)
## Reconstruct
r.ephrinB <- retistruct.reconstruct(o)

plot.flat(r.ephrinB, mesh=FALSE, stitch=FALSE, markup=FALSE)
mtext("G", adj=0, font=2, line=-0.9)
polarplot(r.ephrinB, mesh=FALSE)
mtext("H", adj=0, font=2, line=-0.9)



## To print for PloS Biol.
## dev.copy2eps(file="ChenEtal95.eps", width=6.83, height=6.83/2)
