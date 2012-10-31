\name{trace.ray}
\alias{trace.ray}
\title{Trace ray through system}
\usage{
  trace.ray(r0, S)
}
\arguments{
  \item{r0}{Initial ray vector}

  \item{S}{System defined by data frame containing left
  margin of surfaces (A), radii of surfaces (R), cut-offs
  (C) and refracive indices (n).}
}
\value{
  5-column matrix with Cartesian coordinates of ray origin
  in first two columns (\code{x} and \code{y}), angles of
  rays to x-axis in third column (\code{theta}), index of
  surface at which intersection occurs (\code{j}) and angle
  to x-axis at which intersection occurs on surface
  (\code{alpha}).
}
\description{
  Trace ray through system
}
\author{
  David Sterratt
}

