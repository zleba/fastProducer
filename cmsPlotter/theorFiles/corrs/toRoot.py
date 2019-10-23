#!/usr/bin/env python


import ROOT

import sys
from math import sqrt

def readTable(fName):
    sigmaTab = {}
    #Loop over files
    fp = open(fName, 'r')
    for line in fp:
        line =  line.strip()
        if line == "":
            continue

        l =  map(float, line.split())
        #print l

        y    = l[0]
        ptL  = l[2]
        ptH  = l[3]
        corr = l[4]
        unc = 0

        if "np15" in fName:
            unc = (((l[5] - corr)**2 + (l[6] - corr)**2)/2)**0.5
        elif "np16" in fName:
            unc = ((l[5]**2 + l[6]**2)/2)**0.5 * 0.01 * corr

        if y not in sigmaTab:
            sigmaTab[y] = []
        sigmaTab[y].append( [ptL, ptH, corr, unc] )
        
    return sigmaTab



def rn():
    import random
    return str(random.randint(1,1000000))

def vec(vv):
    from ROOT import std
    import array
    return array.array('d', vv)
    vvv = std.vector("double")()
    for v in vv:
        vvv.push_back(v)
    return vvv



def writeTable(corrTab, name):
    for el in corrTab:
        #get binning
        bins = []
        for x in corrTab[el]:
            bins.append(x[0])
        bins.append(corrTab[el][-1][1])

        hCorr = ROOT.TH1D(rn(), str(el), len(bins)-1, vec(bins))

        #Stat errors
        for i, x in enumerate(corrTab[el]):
            hCorr.SetBinContent(i+1, x[2])
            hCorr.SetBinError(i+1,   x[3])

        hCorr.Write(name+'_y'+str(int(2*el)))
        #print el
        #print bins
        
fOut = ROOT.TFile('np_ew.root', 'RECREATE')

#fileName = "ew15_ak4"

for n in ["ew15_ak4", "ew16_ak4", "ew16_ak7", "np15_ak4", "np16_ak4", "kFactorNLL_ak4", "kFactorNNLO_ak4", "kFactorNLL_ak7", "kFactorNNLO_ak7"]:
    corrTab = readTable(n+'.txt')
    writeTable(corrTab, n)

#print sigmaTab

fOut.Write()
fOut.Close()


import sys
sys.exit()



#print sigmaTab

