##
## Utilities for plotting functions
## 

## Generate colours for strain plots
strain.colours <- function(x) {
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  col <- x/log(0.75)*35 + 35
  col[col<1] <- 1
  col[col>70] <- 70
  return(col)
}

##
## Polar plots
##

## Put text on the polar plot
text.polar <- function(text) {
  mtext(text, 1, adj=1, line=-0.1)
}
