## Installation Requirements (Ubuntu and Mac):

* [R software environment](https://cran.r-project.org/)
* R packages
  * [yaml](https://cran.r-project.org/web/packages/yaml/index.html)
  * [deSolve](http://desolve.r-forge.r-project.org/)
  * [ape](http://ape-package.ird.fr/)
  * [diversitree](https://CRAN.R-project.org/package=diversitree)
  * [RUnit](https://cran.r-project.org/web/packages/RUnit/index.html)
  * [inline](https://cran.r-project.org/web/packages/inline/index.html)
  * [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
  * [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)
  * [whisker](https://cran.r-project.org/web/packages/whisker/index.html)
  * [phyloTop](https://cran.r-project.org/web/packages/phyloTop/index.html) (optional)
  * [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html) (optional)
* GNU tools: 
  * [bison](https://www.gnu.org/software/bison/)
  * [flex](https://github.com/westes/flex)
  * [m4](https://www.gnu.org/software/m4/m4.html) (Only for Mac as it is already installed on Ubuntu)
  * [autoconf](https://www.gnu.org/software/autoconf/autoconf.html) (Only for Mac as it is already installed on Ubuntu)
  * [GFortran](https://gcc.gnu.org/wiki/GFortran) (Only for Mac as it is already installed on Ubuntu)
* C libraries:
  * [GNU Scientific Library](https://www.gnu.org/software/gsl/)
  * [GNU MP Bignum Library](https://gmplib.org/)
  * [Judy C library](http://judy.sourceforge.net/)
  * [C-igraph](http://igraph.org/c/)


## Requirements Installation Procedure (Ubuntu):

* The commands for each step are to be written/copied line by line to the terminal.

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
    install.packages("diversitree")
    install.packages("RUnit")
    install.packages("inline")
    install.packages("Rcpp")
    install.packages("RcppArmadillo")
    install.packages("whisker")
    quit() 
    ```
4. Installing GNU tools and C libraries
    ```
    sudo apt-get install bison
    sudo apt-get install flex  
    sudo apt-get install libgsl-dev
    sudo apt-get install libgmp-dev
    sudo apt-get install libjudy-dev
    sudo apt-get install libigraph0v5
    sudo apt-get install libigraph0-dev
    ```

## Requirements Installation Procedure (Mac):

* The commands for each step are to be written/copied line by line to the terminal.

1. Install the latest version of R from the appropriate [mirror](https://cran.r-project.org/mirrors.html).
2. Installing Xcode
    ```
    xcode-select --install
    ```
   Follow the generated prompts to the end of installation. To verify that Xcode was correctly installed check what version    of Xcode was installed:
    ```
    xcodebuild -version
    ```
3. Installing Command Line Tools

   Go to http://developer.apple.com/downloads and sign in with your Apple ID (the same one you use for iTunes and app
   purchases). Search for "command line tools" (in the search field on the left), then click on version corresponding to the
   installed version of Xcode and click on the the .dmg link to download it. Run the .dmg and follow the generated prompts
   to the end of installation.
4. Installing Homebrew
    ```
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    ```
5. Installing R Packages
    ```
    R
    install.packages("yaml")
    install.packages("deSolve")
    install.packages("ape")
    install.packages("diversitree")
    install.packages("RUnit")
    install.packages("inline") 
    install.packages("Rcpp")
    install.packages("RcppArmadillo")
    install.packages("whisker")
    quit() 
    ```
4. Installing GNU tools and C libraries
    ```
    brew install bison
    brew install flex  
    brew install gsl
    brew install gmp
    brew install m4
    brew install autoconf
    brew install igraph
    ```
    
    GFortran is needed when compiling rcolgem. The .dgm can be download from [here](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS). Run the .dmg and follow the generated prompts to the end of 
    installation.
    
    Judy is no longer supported by Homebrew so we need to compile it from [binaries](https://sourceforge.net/projects/judy/):
    ```
    cd Downloads			# Downloads can be replace by wherever Judy was downloaded.
    tar -xvzf Judy-*.*.*.tar.gz	# The * should be replace by the release number.
    cd judy-*.*.*		# The * should be replace by the release number.
    ./configure
    make
    make check
    sudo make install
    ```
    
## Kaphi Installation Procedure (Ubuntu and Mac):

* Navigate to your preferred location in the filesystem and clone Kaphi from the GitHhub repository
    ```
    git clone --recursive https://github.com/PoonLab/Kaphi.git
    ```
* Compile and install rcolgem
    ```
    cd Kaphi/colgem  # navigate back to the colgem directory
    R CMD INSTALL pkg
    ```
* Compile and install Kaphi
    ```
    cd ..  # navigate back to package root
    R CMD INSTALL pkg
    ```
