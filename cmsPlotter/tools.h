#ifndef tools_H
#define tools_H

#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include <vector>
#include <map>
#include <algorithm>

inline std::vector<double> getRange(double asMin, double asMax, double st = 0.001)
{
    std::vector<double> v;
    for(double as = asMin; as < asMax + st/2; as += st)
        v.push_back(as);
    return v;
}


static const vector<double> ptBinsAs = {97, 174, 272, 395, 548, 737, 967, 1248, 1588, 2000, 2500, 3103};

static const std::map<TString, std::vector<double> > pdfAsVals =  {
    {"CT14nlo",  getRange(0.111, 0.123) },
    {"CT14nnlo", getRange(0.111, 0.123) },

    {"HERAPDF20_NLO",  getRange(0.111, 0.123) },
    {"HERAPDF20_NNLO", getRange(0.111, 0.123) }, 

    {"NNPDF31_nnlo",  {0.112, 0.114, 0.116, 0.117, 0.118, 0.119, 0.120,  0.122}}
};





//Apply NP + EW corrections to theory
inline void applyNPEW(TH1D *h, int y,  TString Year)
{
    TFile *fNPEW  = TFile::Open("theorFiles/corrs/np_ew.root");  //NP+EW corrections

    int tag = Year.Contains("15") ? 15 : 16;

    TH1D *hEW = dynamic_cast<TH1D*>( fNPEW->Get(Form("ew%d_ak4_y%d", tag, y)));
    TH1D *hNP = dynamic_cast<TH1D*>( fNPEW->Get(Form("np%d_ak4_y%d", tag, y)));
    if(!hEW || !hNP) {
        std::cout << "Histogram is missing in np_ew.root file" << std::endl;
        std::exit(0);
    }

    for(int i = 1; i <= h->GetNbinsX(); ++i) {
        double pt  = h->GetBinCenter(i);
        double v   = h->GetBinContent(i);
        double err = h->GetBinError(i);

        int iNP = hNP->FindBin(pt);
        int iEW = hEW->FindBin(pt);

        double np = hNP->GetBinContent(iNP);
        double ew = hEW->GetBinContent(iEW);

        v   *= np * ew;
        err *= np * ew;

        h->SetBinContent(i, v);
        h->SetBinError(i, err);
    }
    fNPEW->Close();
}


//Apply NNLO or NLL k-factor
void applyKfactor(TH1D *h, int y,  TString Tag)
{
    TFile *fNPEW  = TFile::Open("theorFiles/corrs/np_ew.root");  //NP+EW corrections

    TH1D *hCorr = dynamic_cast<TH1D*>( fNPEW->Get(Tag + Form("_ak4_y%d",  y)));
    if(!hCorr) {
        std::cout << "Histogram not found :" << Tag + Form("_ak4_y%d",  y) << std::endl;
        std::exit(0);
    }

    for(int i = 1; i <= h->GetNbinsX(); ++i) {
        double pt  = h->GetBinCenter(i);
        double v   = h->GetBinContent(i);
        double err = h->GetBinError(i);

        int iCorr = hCorr->FindBin(pt);

        double corr = hCorr->GetBinContent(iCorr);

        v   *= corr;
        err *= corr;

        h->SetBinContent(i, v);
        h->SetBinError(i, err);
    }
    fNPEW->Close();
}

TGraphAsymmErrors *getBand(TH1D *hCnt, TH1D *hUp, TH1D *hDn)
{
    hUp->Divide(hCnt);
    hDn->Divide(hCnt);

    TGraphAsymmErrors *gr = new TGraphAsymmErrors(hCnt->GetNbinsX());

    for(int i = 1; i <= hCnt->GetNbinsX(); ++i) {
        double v1 = hUp->GetBinContent(i) - 1;
        double v2 = hDn->GetBinContent(i) - 1;
        double up = std::max(0., std::max(v1, v2));
        double dn = std::max(0., std::max(-v1, -v2));

        gr->SetPoint(i-1, hCnt->GetBinCenter(i), 1);
        gr->SetPointError(i-1, hCnt->GetBinWidth(i)/2., hCnt->GetBinWidth(i)/2., up, dn);
    }
    return gr;
}



#endif
