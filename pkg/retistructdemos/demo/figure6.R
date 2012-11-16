library(retistruct)
## Reconstruct left contra
root <- system.file(package = "retistructdemos")

options(contour.levels=c(25, 50, 75, 95))

if (1) {
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
}

width <- 6.83
height <- width*2/3
## par(mfrow=c(2, 4))
par(mar=c(0.5, 0.5, 1.5, 0.5))
layout(rbind(c(1, 2, 5, 6),
             c(3, 4, 5, 6),
             c(0, 0, 7, 7)),
       widths =c(1, 1, 2, 2),
       heights=c(1, 1, 2))

## 1. Right Ipsi
projection(r.ri, image=FALSE, cex=0.1)
title("Right Ipsi")

## 2. Left Ipsi
projection(r.li, image=FALSE, cex=0.1)
title("Left Ipsi")

## 3. Right Contra
projection(r.rc, image=FALSE, cex=0.1)
title("Right Contra")

## 4. Left Contra
projection(r.lc, image=FALSE, cex=0.1)
title("Left Contra")

## 5. Ipsi orthographic correct
projection(r.li, projection=orthographic,
           transform=invert.sphere.to.hemisphere,
           axisdir=cbind(phi=22, lambda=64),
           proj.centre=cbind(phi=50, lambda=0),
           image=FALSE, pole=TRUE)
projection(r.ri, projection=orthographic,
           transform=invert.sphere.to.hemisphere,
           axisdir=cbind(phi=22, lambda=-64),
           proj.centre=cbind(phi=50, lambda=0),
           image=FALSE, pole=TRUE,
           add=TRUE)
title(expression(paste("Ipsi: O.A. at ", 22 * degree, " el., ", 64 * degree, " az.")))

## 6. Ipsi orthographic incorrect
projection(r.li, projection=orthographic,
           transform=invert.sphere.to.hemisphere,
           axisdir=cbind(phi=35, lambda=60),
           proj.centre=cbind(phi=50, lambda=0),
           image=FALSE, pole=TRUE)
projection(r.ri, projection=orthographic,
           transform=invert.sphere.to.hemisphere,
           axisdir=cbind(phi=35, lambda=-60),
           proj.centre=cbind(phi=50, lambda=0),
           image=FALSE, pole=TRUE,
           add=TRUE)
title(expression(paste("Ipsi: O.A. at ", 35 * degree, " el., ", 60 * degree, " az.")))


## 7. Contra Sinusoidal correct
projection(r.lc, projection=sinusoidal,
           transform=invert.sphere.to.hemisphere,
           axisdir=cbind(phi=22, lambda=64),
           proj.centre=cbind(phi=0, lambda=0),
           cex=0.1,
           image=FALSE, pole=TRUE)
projection(r.rc, projection=sinusoidal,
           transform=invert.sphere.to.hemisphere,
           axisdir=cbind(phi=22, lambda=-64),
           proj.centre=cbind(phi=0, lambda=0),
           cex=0.1,
           image=FALSE, pole=TRUE,
           add=TRUE)
title(expression(paste("Contra: O.A. at ", 22 * degree, " el., ", 64 * degree, " az.")))


dev.print(pdf, "figure6.pdf", width=width, height=height)


