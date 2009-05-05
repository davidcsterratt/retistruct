## Copyright (C) 2007 David Bateman
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{idx}, @var{p}] =} tsearchn (@var{x}, @var{t}, @var{xi})
## Searches for the enclosing Delaunay convex hull. For @code{@var{t} =
## delaunayn (@var{x})}, finds the index in @var{t} containing the
## points @var{xi}. For points outside the convex hull, @var{idx} is NaN.
## If requested @code{tsearchn} also returns the barycentric coordinates @var{p}
## of the enclosing triangles.
## @seealso{delaunay, delaunayn}
## @end deftypefn

tsearchn <- function(x, t, xi) {
##  if (nargin != 3)
##    print_usage ();
##  endif

  nt <- dim(t)[1]
  m <- dim(x)[1]
  n <- dim(x)[2]
  mi <- dim(xi)[1]
  idx <- rep(NA, mi)
  p <- matrix(NA, mi, n + 1)

  ni <- 1:mi
  for (i in 1:nt) { 
    ## Only calculate the Barycentric coordinates for points that have not
    ## already been found in a triangle.
    b <- cart2bary(x[t[i,],], xi[ni,,drop=FALSE]);

    ## Our points xi are in the current triangle if
    ## (all(b >= 0) && all (b <= 1)). However as we impose that 
    ## sum(b,2) == 1 we only need to test all(b>=0). Note need to add
    ## a small margin for rounding errors
    intri <- apply(b >= -1e-12, 1, all)
    if ((sum(intri))!=0) {
      print(paste(sum(intri), "points found in triangle", i))
    }
    idx[ni[intri]] <- i
    p[ni[intri],] <- b[intri, ]
    ## ni[intri] <- FALSE;
    ni <- ni[!intri]
  }
  return(list(idx=idx, p=p))
}

cart2bary <- function(T, P) {
  ## Conversion of Cartesian to Barycentric coordinates.
  ## Given a reference simplex in N dimensions represented by a
  ## (N+1)-by-(N) matrix, and arbitrary point P in cartesion coordinates,
  ## represented by a N-by-1 row vector can be written as
  ##
  ## P = Beta * T
  ##
  ## Where Beta is a N+1 vector of the barycentric coordinates. A criteria
  ## on Beta is that
  ##
  ## sum (Beta) == 1
  ##
  ## and therefore we can write the above as
  ##
  ## P - T(end, :) = Beta(1:end-1) * (T(1:end-1,:) - ones(N,1) * T(end,:))
  ##
  ## and then we can solve for Beta as
  ##
  ## Beta(1:end-1) = (P - T(end,:)) / (T(1:end-1,:) - ones(N,1) * T(end,:))
  ## Beta(end) = sum(Beta)
  ##
  ## Note below is generalize for multiple values of P, one per row.
  M <- dim(P)[1]
  N <- dim(P)[2]
  ## 
  Beta <- (P - matrix(T[N+1,], M, N, byrow=TRUE)) %*% solve(T[1:N,] - matrix(1,N,1) %*% T[N+1,,drop=FALSE])
  Beta <- cbind(Beta, 1 - apply(Beta, 1, sum))
  return(Beta)
}

bary2cart <- function(T, Beta) {
  ## Conversion of Barycentric to Cartesian coordinates.
  ## Given a reference simplex T in N dimensions represented by a
  ## (N+1)-by-(N) matrix, and arbitrary point Beta in baryocentric coordinates,
  ## represented by a N+1-by-1 row vector, the cartesian coordinates P are
  ## given
  ##
  ## P = Beta * T
  return(Beta %*% T)
}
