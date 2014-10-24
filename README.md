Installation
============

* To install the files provided with the PLoS Comp. Biol. paper under
  GNU/Linux and Mac OS X, look at the the user guide
  (doc/retistruct-user-guide.pdf) for how to install R and
  dependencies, but do the installation from within R by typing:

    source("install.R")

  and

    source("install-gui.R")

  instead of the source statements in the manual. There are no binary
  packages available for Windows; to build Retistruct from source,
  follow the instructions at
  http://www.biostat.wisc.edu/~kbroman/Rintro/Rwinpack.html

* To install the latest version and for user instructions, see
  doc/retistruct-user-guide.pdf

* For description of functions within the package, see
  doc/retistruct-manual.pdf

Troubleshooting
===============

If you get an error like this:

  Error in asCairoDevice() : Graphics API version mismatch

in Ubuntu, then try the following at the command line:

1) At the bash prompt:
sudo apt-get install build-essential libgtk2.0-dev

2) At the R prompt:
install.packages("cairoDevice")

(See https://bugs.launchpad.net/ubuntu/+source/cairodevice/+bug/778804)

