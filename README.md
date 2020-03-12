# <code>ipbridging</code> Package

## Authors

Tzu-Ping Liu (UC Davis) <br>
Gento Kato (UC Davis) <br>
Sam Fuller (UC Davis)

## Description

Bridging Ideal Point Estimates. This package provides parametric and non-parametric methods to bridge ideal point estimates generated from two (or more) separate data sets. 
The package website is published [HERE](https://gentok.github.io/ipbridging/).

## Installation

### 1. Insall <code>oc</code> package

<b>CAUTION!</b> As of January 2020, one of the dependencies, <code>oc</code> package is archived from CRAN. 
If you don't already have <code>oc</code> package, you can install the latest version of archived <code>oc</code> package using following codes:

<code>install.packages("https://cran.r-project.org/src/contrib/Archive/oc/oc_1.01.tar.gz", repos=NULL, type="source")</code>

Alternatively, you could install <code>oc</code> package from github CRAN mirror:

<code>devtools::install_github("cran/oc")</code>

### 2. Install <code>ooc</code> package

<b>CAUTION!</b> May need the first line to avoid error in installation due to warnings.

<code>Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)</code> <br>
<code>devtools::install_github("tzuliu/ooc")</code>

### 3. Install <code>ipbridging</code> package

After successfully installing <code>oc</code> and <code>ooc</code> package, you can install the <code>ipbridging</code> package by:

<code>Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)</code> <br>
<code>devtools::install_github("gentok/ipbridging")</code>

<!-- Also, in case you want to use GPU accelerated method to tune SVM parameters in <code>oocflex</code> function, you need to install <code>Rgtsvm</code> package (check the GitHub repository [HERE](https://github.com/Danko-Lab/Rgtsvm)). You need to have PC with NVIDIA GPU and Linux OS to install <code>Rgtsvm</code> package. -->

## Main Functions

* <code>[ipbridging](https://gentok.github.io/ipbridging/reference/ipbridging.html)</code>: The main function from this package. It implements methods to bridge ideal point estimates generated from two data sets.
* <code>[oocflex](https://gentok.github.io/ipbridging/reference/oocflex.html)</code>: The implementation of Ordered Optimal Classification (OOC) with flexible estimation strategies. This function is the modified version of <code>ooc</code> function in [ooc package](https://github.com/tzuliu/ooc).

## Updates Log

* 03/12/2020 Version 0.0.08 (beta version 8: Made Significant Changes)
* 02/03/2020 Version 0.0.06 (beta version 7)
* 01/29/2020 Version 0.0.05 (beta version 6)
* 01/29/2020 Version 0.0.04 (beta version 5)
* 01/22/2020 Version 0.0.03 (beta version 4)
* 01/21/2020 Version 0.0.02 (beta version 3)
* 01/20/2020 Version 0.0.01 (beta version 2)
* 01/19/2020 Version 0.0.00 (beta version) released
