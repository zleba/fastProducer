#!/usr/bin/env python

import ROOT

import sys
from math import sqrt

def readTable():
    sigmaTab = {}
    #Loop over files
    fp = open('13TeVorg.txt', 'r')
    n = 0
    arr = []
    isIn = False
    for line in fp:
        line =  line.strip()
        if line == "*":
            isIn = True
            continue
        if not isIn:
            continue

        l =  map(float, line.split())
        #print l

        y    = l[1]
        ptL  = l[3]
        ptH  = l[4]
        sigma= l[5]
        stat = l[6]
        uncor= l[7]

        sys= l[10:]
        sysH = sqrt(sum([max(0, x)**2 for x in sys]))
        sysL = sqrt(sum([max(0,-x)**2 for x in sys]))

        if y not in sigmaTab:
            sigmaTab[y] = []
        sigmaTab[y].append( [ptL, ptH, sigma, 0.01 * (stat**2+uncor**2)**0.5 * sigma, 0.01*sysH*sigma, 0.01*sysL*sigma] )
        
    return sigmaTab

sigmaTab = readTable()

print sigmaTab

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


fOut = ROOT.TFile('dataEarly.root', 'RECREATE')

for el in sigmaTab:
    #get binning
    bins = []
    for x in sigmaTab[el]:
        bins.append(x[0])
    bins.append(sigmaTab[el][-1][1])

    hStat = ROOT.TH1D(rn(), str(el), len(bins)-1, vec(bins))
    hSysUp = hStat.Clone(rn())
    hSysDn = hStat.Clone(rn())

    #Stat errors
    for i, x in enumerate(sigmaTab[el]):
        hStat.SetBinContent(i, x[2])
        hStat.SetBinError(i,   x[3])

    for i, x in enumerate(sigmaTab[el]):
        hSysUp.SetBinContent(i, x[2] + x[4])
        hSysDn.SetBinContent(i, x[2] - x[4])

    hStat.Write('hStat_y'+str(int(2*el)))
    hSysUp.Write('hSysUp_y'+str(int(2*el)))
    hSysDn.Write('hSysDn_y'+str(int(2*el)))
    print el
    print bins
        
fOut.Write()
fOut.Close()
#print sigmaTab

