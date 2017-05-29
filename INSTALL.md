## Installation Requirements:

* [R software environment](https://cran.r-project.org/)
* R packages
  * [yaml](https://cran.r-project.org/web/packages/yaml/index.html)
  * [deSolve](http://desolve.r-forge.r-project.org/)
  * [ape](http://ape-package.ird.fr/)
  * [diversitree](https://CRAN.R-project.org/package=diversitree)
  * [ECcmtc](https://cran.r-project.org/web/packages/ECctmc/index.html)
  * [RUnit](https://cran.r-project.org/web/packages/RUnit/index.html)
  * [inline](https://cran.r-project.org/web/packages/inline/index.html)
  * [rcolgem](http://colgem.r-forge.r-project.org/)
* GNU tools: 
  * [bison](https://www.gnu.org/software/bison/)
  * [flex](https://github.com/westes/flex)
* C libraries:
  * [GNU Scientific Library](https://www.gnu.org/software/gsl/) 
  * [Judy C library](http://judy.sourceforge.net/) 


## Requirements Installation Procedure (Ubuntu):

* The commands for each step are to be written/coppied one by one to the terminal.

1. Updating and Upgrading The System  
    ```
       sudo apt-get update
       sudo apt-get upgrade
    ```
2. Installing R
    ```
	   sudo apt-get install r-base
	   sudo apt-get install r-base-dev
    ```
3. Installing R Packages
    ```
	   R
	   install.packages("yaml")
	   install.packages("deSolve")
	   install.packages("ape")
	   install.packages("ECcmtc")
	   install.packages("RUnit")
	   install.packages("inline")
	   install.packages("rcolgem", repos="http://R-Forge.R-project.org") 
	   quit() 
    ```
4. Installing GNU tools and C libraries
    ```
	   sudo apt-get install bison
	   sudo apt-get install flex  
	   sudo apt-get install libgsl-dev
	   sudo apt-get install libjudy-dev
    ```
    
## Kaphi Installation Procedure (Ubuntu):

* Navigate to your preferred location in the filesystem and clone Kaphi from the GitHhub repository
    ```
	   git clone --recursive https://github.com/PoonLab/Kaphi.git
    ```
    
* Compile igraph
    ```
	   cd Kaphi/pkg/src/igraph
	   touch configure.ac aclocal.m4 configure Makefile.am Makefile.in
	   ./configure
	   make
	   sudo make install
    ```
* Compile and install Kaphi
    ```
	   cd ../../../  # navigate back to package root
	   R CMD INSTALL pkg
    ```
