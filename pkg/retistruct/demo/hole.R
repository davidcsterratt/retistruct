oldpar <- par(no.readonly=TRUE) # Save graphics parameters before plotting

par(mfcol=c(1, 2))

P <- list(rbind(c(1,1.5),
                c(1.5,1),
                c(2,1),
                c(2.5,2),
                c(3,1),
                c(4,1),
                c(1,4)),
          rbind(c(-1.5,1),
                c(-1,1.5),
                c(-1,4),
                c(-2,3),
                c(-2,2),
                c(-3,2),
                c(-4,1)),
          rbind(c(-4,-1),
                c(-1.5,-1),
                c(-1,-1.5),
                c(-1,-4)),
          rbind(c(1,-1.5),
                c(1.5,-1),
                c(2,-1),
                c(2.5,-2),
                c(3,-1),
                c(4,-1),
                c(1,-4)))
## Stitched outlines
a <- StitchedOutline$new(P)

## Set a fixed point
## One that is in the rim should be fine
a$setFixedPoint(6, "Nasal")

## One that is not in the rim should be moved
a$addTear(c(11, 12, 13))
a$addTear(c(3, 4, 5))
a$addTear(c(21, 22, 23))

## Add fullcuts
a$addFullCut(c(2, 6, 20, 24))
a$addFullCut(c(1, 7, 9, 10))
a$addFullCut(c(8, 14, 15, 16))
a$addFullCut(c(17, 18, 19, 25))

r <- ReconstructedOutline$new()
r$loadOutline(a)

flatplot(a)
flatplot(r)

r$reconstruct()

flatplot(r)
projection(r)

par(oldpar) # Restore graphics parameters
