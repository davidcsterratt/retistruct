* For the main installation and user instructions, see
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
