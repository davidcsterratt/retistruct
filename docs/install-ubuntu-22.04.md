# Howto install Retistruct on Ubuntu 22.04

## 0. Get Ubuntu 22.04

Either:
1. Install Ubuntu on your machine as per the instructions at
   https://ubuntu.com/desktop
2. Install a VirtualBox virtual machine and install and run Ubuntu
   22.04 within this machine:
   https://linux.how2shout.com/how-to-install-ubuntu-22-04-lts-iso-in-virtualbox-vm-to-test-it/

## 1. Install the latest version of R and the required libraries

1. Log into your Ubuntu machine
2. Open a Terminal window
3. Type the following commands

```
# R and relevant R packages
# Taken from https://cran.r-project.org/bin/linux/ubuntu/
sudo apt install --no-install-recommends software-properties-common dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt install r-base r-base-dev
sudo add-apt-repository ppa:c2d4u.team/c2d4u4.0+
sudo apt install --no-install-recommends r-cran-devtools

# Other packages required to install Retistruct
sudo apt-get install libgtk2.0-dev libgl1-mesa-dev libglu1-mesa-dev
```

## 2. Install the development version of Retistruct

1. Open a Terminal
2. Type `R` and Return
3. Paste in
```
devtools::install_github('https://github.com/lawremi/cairoDevice')
devtools::install_github('https://github.com/lawremi/RGtk2/RGtk2', subdir='RGtk2')
devtools::install_github('https://github.com/jverzani/gWidgets2RGtk2', force=TRUE)
devtools::install_github("davidcsterratt/retistruct@v0.7.x", subdir="pkg/retistruct")
```
4. Now try typing:
```
library(retistruct)
retistruct()
```
