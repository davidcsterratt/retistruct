library(retistruct)
## Reconstruct left contra
root <- system.file(package = "retistructdemos")

dataset <- file.path(root, "extdata", "Figure_6-data", "left-contra")
r <- retistruct.read.dataset(dataset)
r <- retistruct.read.markup(r)
r.lc <- retistruct.reconstruct(r)

## Reconstruct right ipsi
dataset <- file.path(root, "extdata", "Figure_6-data", "left-ipsi")
r <- retistruct.read.dataset(dataset)
r <- retistruct.read.markup(r)
r.li <- retistruct.reconstruct(r)

## Reconstruct right contra
dataset <- file.path(root, "extdata", "Figure_6-data", "right-contra")
r <- retistruct.read.dataset(dataset)
r <- retistruct.read.markup(r)
r.rc <- retistruct.reconstruct(r)

## Reconstruct right ipsi
dataset <- file.path(root, "extdata", "Figure_6-data", "right-ipsi")
r <- retistruct.read.dataset(dataset)
r <- retistruct.read.markup(r)
r.ri <- retistruct.reconstruct(r)

width <- 6.83
height <- width/2.5
## par(mfrow=c(2, 4))
par(mar=c(0.5, 0.5, 1.5, 0.5))
layout(matrix(1:8, 2, 4, byrow=TRUE), widths=c(1, 2, 2, 1))
## Right Contra
projection(r.rc)
projection(r.rc, projection=sinusoidal, transform=invert.sphere, axisdir=cbind(phi=145, lambda=30))
title("Right Contra")

## Left Contra
projection(r.lc, projection=sinusoidal, transform=invert.sphere, axisdir=cbind(phi=145, lambda=150))
title("Left Contra")
projection(r.lc)


## Right Ipsi
projection(r.ri)
projection(r.ri, projection=sinusoidal, transform=invert.sphere, axisdir=cbind(phi=145, lambda=30))
title("Right Ipsi")


## Left Ipsi
projection(r.li, projection=sinusoidal, transform=invert.sphere, axisdir=cbind(phi=145, lambda=150))
title("Left Ipsi")
projection(r.li)

dev.print(pdf, "figure6.pdf", width=width, height=height)
