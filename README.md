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

To install the current stable verision (i.e. the one currently on CRAN), use the R devtools package like this:
```
devtools::install_github("davidcsterratt/retistruct@v0.6.1", subdir="pkg/retistruct")
```
To install from the stable development branch (i.e. code that will be in the next CRAN release), use the R devtools package like this:
```
devtools::install_github("davidcsterratt/retistruct", subdir="pkg/retistruct")
```

Roadmap
=======

There are a number of improvements on the horizon - see the [list of milestones](https://github.com/davidcsterratt/retistruct/milestones) for full details.

Retistruct is not my main work at present, so improvements may take some time. Nevertheless, if there's something you'd like implemented, please let me know, ideally by [creating an issue](https://github.com/davidcsterratt/retistruct/issues/new).

Funding Acknowldegements
========================

The development of the initial version of Retistruct was supported by
a Programme Grant from the UK Wellcome Trust (G083305) from 2008-2013.

Improvements to image handing and refactoring the code (in the
"master" github branch and to be released in v0.6.0) were supported by
The Jackson Laboratory (Bar Harbor, ME, USA) Scientific Services
Innovation Fund from 2016-2017 and an NIH R21 grant (EY027894–01A1)
from 2018-2020 to Dr. Mark P. Krebs, The Jackson Laboratory.

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/davidcsterratt/retistruct.svg?branch=master)](https://travis-ci.org/davidcsterratt/retistruct)
<!-- badges: end -->

