## Installing R on Ubuntu
Installing R on Ubuntu may require a couple of extra steps compared to Windows or Mac OS (notably to allow non-base R packages to be installed). These notes document one way to get R up and running. The steps were performed on a clean install of Ubuntu v14.04 running through a VirtualBox (v5.0.12 r104815) Virtual Machine on Mac OS 10.11. I suspect there might be more efficient ways to do this, but wanted to document one way that worked for me. 

These steps have not been tested on other Linux systems. Further instructions are provided under the "Download R for Linux" link and README files for each Linux OS <a href="https://cran.r-project.org/" target="_blank">available through the CRAN website</a>. This is a distilled version of those instructions for Ubuntu, with a couple added commands that I found useful in my experience. 

-------

Start by editing the `/etc/apt/sources.list` file. You can use the `gedit` text editor that comes pre-installed on Ubuntu, but you will need administrative privileges to save changes to the `sources.list` file. Open the terminal and type this command to open `gedit` with administrative privileges.

```shell
sudo gedit /etc/apt/sources.list
```
Add this line to the end of the file (maybe with a '# Comment noting this is for R'), save and close it.

    deb http://cran.rstudio.com/bin/linux/ubuntu trusty/

Back at the terminal, add the public key for Ubuntu on CRAN, which will help when installing R packages:
```shell
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
```

Now we can update the available list of packages. This will ensure that the latest version of R is installed during the next step (otherwise, Ubuntu may install an older version).
```shell
sudo apt-get update
```

To install R (this could take a couple of minutes), open the terminal and type:
```shell
sudo apt-get -y install r-base 
```

I recommend running R through RStudio, a free IDE, available at: https://www.rstudio.com/

Also, I recommend starting RStudio (or just R) from the terminal using `sudo`. This will ensure that the R session has admin privileges to install packages.
```shell
sudo rstudio
# Alternatively:
# sudo R
```

