# Jet data vs NLO plotter

## General setup

To get the plots with predictions first setup environment from cvmfs
```
. ./setup.sh
```
next install fastNLO and plottingHelper:
```
./installFast.sh
./installPlHelper.sh
```

## Getting theory

Go to directory `cmsPlotter`.
Compile the code.
```
make calcTheory
```
Run the fastNLO
```
./calcTheory
```

## Converting data from xFitter-text files to root files


## Plotting histograms
