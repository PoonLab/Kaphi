## Installation Requirements:

* R Programming language
* R Pakages: yaml, deSolve, ape, ECcmtc, RUnit
* Linux Tools: bison, flex, libgsl-dev, libjudy-dev


## Requirements Installation Procedure (Ubuntu):

* The commands for each step are to be written/coppied one by one to the terminal.

* Updating and Upgrading The System  
    ```
       cd
       sudo apt-get update
       sudo apt-get upgrade
    ```
* Installing R
    ```
	   cd
	   sudo apt-get install r-base
    ```
* Installing R Packages
    ```
	   cd
	   R
    ```
    * In CRAN run the following:    
    ```
	   install.packages("yaml")
	   install.packages("deSolve")
	   install.packages("ape")
	   install.packages("ECcmtc")
	   install.packages("RUnit")
    ```
    * Exit CRAN by running:
    ```
	   quit() 
    ```
* Installing Linux Tools
    ```
	   cd
	   sudo apt-get install bison
	   sudo apt-get install flex  
	   sudo apt-get install libgsl-dev
	   sudo apt-get install libjudy-dev
    ```
    
## Kaphi Installation Procedure (Ubuntu):

* Clone Kaphi In Your Desired Git Directory (for this example said directory will be home)
    ```
	   cd	   
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
* Install Kaphi
    ```
	   cd ~/Kaphi
	   R CMD INSTALL pkg
    ```
