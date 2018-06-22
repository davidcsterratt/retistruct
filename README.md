Retistruct
==========

_computational reconstruction and transformation of flattened retinae_

Retistruct is an R package to morph a flat surface with cuts (a
dissected flat-mount retina) onto a curvilinear surface (the a
standard retinal shape). It can estimate the position of a point on
the intact adult retina to within 8° of arc (3.6% of nasotemporal
axis). The coordinates in reconstructed retinae can be transformed to
visuotopic coordinates.

For full details go to the home page: http://davidcsterratt.github.io/retistruct/

Installation from Github
========================

Most users will wish to install from CRAN, as described on the [Retistruct homepage](http://davidcsterratt.github.io/retistruct/).

To install from the current stable branch (i.e. the branch that is on CRAN), use the R devtools package like this:
```
devtools::install_github("davidcsterratt/retistruct@v0.5.x", subdir="pkg/retistruct")
```
To install from the current, unstable, lesser tested development branch, use the R devtools package like this:
```
devtools::install_github("davidcsterratt/retistruct", subdir="pkg/retistruct")
```

Roadmap
=======

Here are a number of improvements on the horizon. Retistruct is not my main work at present, so improvements may take some time. Nevertheless, if there's something you'd like implemented, please let me know.

* v0.6.x: Major refactor using R6 classes, much improved stitching algorithm
* v0.7.x: Automatic rim angle determination

