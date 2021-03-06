# Jet data vs NLO plotter & Fitter

## General setup

To get the plots with predictions first setup environment from cvmfs
```
. ./setup.sh
```
next install fastNLO and plottingHelper (just for the first time):
```
./installFast.sh
./installPlHelper.sh
```

## Getting theory

The general idea is to store the theory predictions calculated by the fastNLO (from table files cmsPlotter/theorFiles/\*.tab) to the root file histograms.
This can be done as.

Go to directory `cmsPlotter`.
Compile the code.
```
make calcTheory
```
Run the fastNLO
```
./calcTheory
```

We keep two versions of the theory:
1) For the comparison with data, includes:
- scale unc. band
- PDF unc. band
- alphaS unc. band (+- 0.0015 according to PDF4LHC paper)
In all cases the unc is stored simply by two histograms - one for up and one down variation

2) For the alphaS fitting, which includes many theory variants
- for all alphaS values where PDFs are available
- for each PDF value all 7 scale variations are provided
- in case of default alphaS (0.118) also predictions for all PDF eigenvectors are stored

In principle one can derive 1) from 2).
The theory does not include NP&EW corrections and possible k-factors, since these can be fastely applied before plotting or fitting.

All this is provided for several PDFs at NLO & NNLO and all the rapidity bins of the analysis.


## Converting data from xFitter-text files to root files
For plotting the xFitter files can be converted to the root files, this can be done using
```
cmsPlotter/xFitterTables/toRoot.py
```
The converted version is suitable only for plotting, since it contains only stat & overall systematic unc.


## Plotting histograms
For plotting use the macro
```
cmsPlotter/plotJets.C
```


## Fitting histograms
This is about fits of the alphaS only which are based on some global PDF.
One can also calculate the chi2 of the nominal theory (with aS=0.118) to get the chi2 for the given PDF
In addtion to extraction of the aS value, the smoothnes of the data can be tested by a plotting tool, to check theory/data consistency.

In the fits, all unc. are considered, including theory unc. from PDFs nuisance parameters and unc. of data (sys shifts + stat errors)
The main code is in file
```
cmsPlotter/fitTheory.cc
```

