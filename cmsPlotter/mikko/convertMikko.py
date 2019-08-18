#!/usr/bin/env python

import ROOT

def rn():
    import random
    return str(random.randint(1,101000))

#second is template
def rebin(h, hTemp):
    hNew =  hTemp.Clone(rn())
    hNew.Reset()

    for i in range(1, h.GetNbinsX() + 1):
        xCnt = h.GetBinCenter(i)
        iBin = hNew.FindBin(xCnt)

        if(iBin == 0 or iBin == hNew.GetNbinsX()+1):
            continue

        v  = h.GetBinContent(i)
        er = h.GetBinError(i)
        w  = h.GetBinWidth(i)

        orgVal = hNew.GetBinContent(iBin)
        orgErr = hNew.GetBinError(iBin)
        hNew.SetBinContent(iBin, orgVal + v*w)
        from math import hypot
        hNew.SetBinError(iBin, hypot(orgErr,  er*w))


    hNew.Scale(1, "width")
    return hNew





fDet  = ROOT.TFile("common2016_V11_hotzone-3.root")
fCorr = ROOT.TFile("unfold.root")
fMet  = ROOT.TFile("Pythia16Flat_forMikko.root")

yBins = ["0.0-0.5", "0.5-1.0", "1.0-1.5", "1.5-2.0", "2.0-2.5"]

hPart = []
for y in range(5):
    hDet =  fDet.Get("ak4/Eta_" + yBins[y] + "/hpt_data_2016_all_det")
    hCorr =  fCorr.Get("hr_" + yBins[y] + "_2016")
    hP = hDet.Clone("hStat_y" + str(y))
    hP.Divide(hCorr)
    hPart.append(hP)

fTemp = ROOT.TFile.Open("../xFitterTables/data16.root")
hTemp = [fTemp.Get("hStat_y" + str(i)) for i in range(5)] 

for i in range(5):
    hPart[i] = rebin(hPart[i], hTemp[i])
    hBef = fMet.Get("before").ProjectionX(rn(), i+1, i+1)
    hAft = fMet.Get("after").ProjectionX(rn(), i+1, i+1)
    hBef = rebin(hBef, hTemp[i])
    hAft = rebin(hAft, hTemp[i])
    hAft.Divide(hBef)

    if i == 0:
        hAft.Print("all")

    hPart[i].Divide(hAft)

#for i in range(1, hPart[0].GetNbinsX()+1):
#    print "Mikko ", hPart[0].GetBinLowEdge(i), hPart[0].GetBinContent(i)
#
#for i in range(1, hTemp.GetNbinsX()+1):
#    print "Ours ", hTemp.GetBinLowEdge(i), hTemp.GetBinContent(i)


fOut = ROOT.TFile.Open("../xFitterTables/data16m.root", "RECREATE")
for i in range(5):
    hPart[i].Write("hStat_y" + str(i))
fOut.Write()
fOut.Close()
