Retistruct
==========

_computational reconstruction and transformation of flattened retinae_

Retistruct is an R package to morph a flat surface with cuts (a
dissected flat-mount retina) onto a curvilinear surface (a
standard retinal shape). It can estimate the position of a point on
the intact adult retina to within 8° of arc (3.6% of nasotemporal
axis). The coordinates in reconstructed retinae can be transformed to
visuotopic coordinates.

For full details go to the home page: http://davidcsterratt.github.io/retistruct/

Installation from Github
========================

Retistruct is currently not installable from CRAN, as described on the
[Retistruct homepage](http://davidcsterratt.github.io/retistruct/). In
addition, a number of the packages on which Retistruct depends are no
longer on CRAN. Therefore, you need to install from Github.

## Install the `devtoools` package to enable installation from Github

First, from the R console, install the `devtools` package by typing:
```
install.packages("devtools")
```
followed by `Return`. This command will take some time to run, as it has to install a number of packages that it depends on.

## Install the requirements for the graphical user interface (GUI)

To use the graphical user interface of Retistruct, you will need to
install the GTK libraries as follows:

* Mac: using [this guidance](https://github.com/davidcsterratt/retistruct/issues/4)
* Ubuntu Linux:
  ```
  sudo apt-get install r-base r-base-dev libgtk2.0-dev libgl1-mesa-dev libglu1-mesa-dev
  ```

Then install three R packages, as follows:
```
devtools::install_github('https://github.com/lawremi/cairoDevice')
devtools::install_github('https://github.com/lawremi/RGtk2/RGtk2', subdir='RGtk2')
devtools::install_github('https://github.com/jverzani/gWidgets2RGtk2', force=TRUE)
```

## Install Retistruct

Then, do one of the following from the R console:

* To install the latest development branch (not yet released on
   CRAN):
   ```
   devtools::install_github("davidcsterratt/retistruct@v0.7.x", subdir="pkg/retistruct")
   ```
   The development [User guide](https://github.com/davidcsterratt/retistruct/blob/v0.7.x/docs/retistruct-user-guide.pdf) has more details on how to use the new features.

* To install from the stable development branch, use the R devtools package like this:
   ```
   devtools::install_github("davidcsterratt/retistruct", subdir="pkg/retistruct")
   ```

* To install the current stable version, use the R devtools package
   like this: ```
   devtools::install_github("davidcsterratt/retistruct@v0.6.4",
   subdir="pkg/retistruct") ``` You will need to replace `0.6.4` with
   the version number of the latest release.

Roadmap
=======

There are a number of improvements on the horizon - see the [list of milestones](https://github.com/davidcsterratt/retistruct/milestones) for full details.

Retistruct is not my main work at present, so improvements may take some time. Nevertheless, if there's something you'd like implemented, please let me know, either by [creating an issue](https://github.com/davidcsterratt/retistruct/issues/new), or [by email](mailto:david.c.sterratt@ed.ac.uk).

Funding Acknowledgements
========================

The development of the initial version of Retistruct was supported by
a Programme Grant from the UK Wellcome Trust (G083305) from 2008-2013.

Improvements to image handing and refactoring the code (released in
v0.6.0) were supported by The Jackson Laboratory (Bar Harbor, ME, USA)
Scientific Services Innovation Fund from 2016-2017 and an NIH R21
grant (EY027894) from 2018-2020 to Dr. Mark P. Krebs, The Jackson
Laboratory.

The capabilities to reconstruct tissue comprised of separate fragments
(released in v0.7.0) and to reconstruct 3D data comprising an overhead
image and depth map (released in v0.7.2), and user interface
improvements (released in v0.7.0) were supported by an NIH R21 grant
(EY027894) from 2018-2020 to Dr. Mark P. Krebs, The Jackson
Laboratory.

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/davidcsterratt/retistruct.svg?branch=master)](https://travis-ci.com/davidcsterratt/retistruct)
<!-- badges: end -->


<!--  LocalWords:  Retistruct Github CRAN devtools davidcsterratt EY
 -->
<!--  LocalWords:  subdir retistruct Roadmap Wellcome Harbor
 -->
