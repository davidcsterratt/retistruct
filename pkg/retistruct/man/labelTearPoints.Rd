\name{labelTearPoints}
\alias{labelTearPoints}
\title{Label three outline point indicies as apicies and vertices of tear}
\usage{
labelTearPoints(o, m)
}
\arguments{
  \item{m}{the vector of three indicies}

  \item{o}{Outline object}
}
\value{
Vector of indicies labelled with V0, VF and VB
}
\description{
Label a set of three unlabelled points supposed to refer to
the apex and vertcies of a cut and tear with the V0 (Apex),
VF (forward vertex) and VB (backward vertex) labels.
}
\author{
David Sterratt
}

