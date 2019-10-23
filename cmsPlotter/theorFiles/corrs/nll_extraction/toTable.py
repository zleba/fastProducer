#!/usr/bin/env python

def getVals(fName):
    f = open(fName, "r")

    sigmaTab = {}

    for x in f:
        x = x.strip()
        if len(x) == 0: continue
        if x[0].isdigit() == False: continue
        x = [x for x in x.split(' ') if x != '']
        #print(x)
        #print(x[0], x[2])
        #y == 0

        for a in (4, 9, 14, 19, 24):
            if a+1 >= len(x): continue
            if '-----' in x[a]: continue
            m = x[a].replace('$\\,\\cdot\\,', '')
            e = x[a+1].replace('10^{', '').replace('}$', '')
            #print ( m + 'e' +  e)

            y = (a-4)/5
            ptL = int(x[0])
            ptH = int(x[2])
            sigma = float(m+'e'+e)
            if y not in sigmaTab:
                sigmaTab[y] = []
            sigmaTab[y].append( [ptL, ptH, sigma] )
    return sigmaTab

def divideTabs(sigmaTab, sigmaTabNll):
    assert(len(sigmaTab) == len(sigmaTabNll))

    corrTab = {}

    for y in sigmaTab:
        assert(len(sigmaTab[y]) == len(sigmaTabNll[y]))
        corrTab[y] = []

        for i in range(len(sigmaTab[y])):
            assert(sigmaTab[y][i][0] == sigmaTabNll[y][i][0] )
            assert(sigmaTab[y][i][1] == sigmaTabNll[y][i][1] )
            ptL = sigmaTab[y][i][0]
            ptH = sigmaTab[y][i][1]
            c1  = sigmaTab[y][i][2]
            c2  = sigmaTabNll[y][i][2]

            corrTab[y].append( [ptL, ptH, c2 / c1] )
    return corrTab

def printToFile(corrTab):
    for y in corrTab:
        for l in corrTab[y]:
            print y*0.5, (y+1)*0.5, l[0], l[1], l[2]

R = 7
sigmaTab    = getVals("table-cms13TeV-R0"+str(R)+".tex")
sigmaTabNll = getVals("table-cms13TeV-R0"+str(R)+"-nll.tex")

corr = divideTabs(sigmaTab, sigmaTabNll)
#print corr
printToFile(corr)
