R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)

#include "plottingHelper.h"
using namespace PlottingHelper;//pollute the namespace!

TString rn() {return Form("%d",rand());}

TString year = "15";


TGraphAsymmErrors *getBand(TH1D *hCnt, TH1D *hUp, TH1D *hDn)
{
    hUp->Divide(hCnt);
    hDn->Divide(hCnt);

    TGraphAsymmErrors *gr = new TGraphAsymmErrors(hCnt->GetNbinsX());

    for(int i = 1; i <= hCnt->GetNbinsX(); ++i) {
        double v1 = hUp->GetBinContent(i) - 1;
        double v2 = hDn->GetBinContent(i) - 1;
        double up = max(0., max(v1, v2));
        double dn = max(0., max(-v1, -v2));

        gr->SetPoint(i-1, hCnt->GetBinCenter(i), 1);
        gr->SetPointError(i-1, hCnt->GetBinWidth(i)/2., hCnt->GetBinWidth(i)/2., up, dn);
    }
    return gr;
}

//Apply NP + EW corrections to theory
void applyNPEW(TH1D *h, int y,  TString Year)
{
    TFile *fNPEW  = TFile::Open("theorFiles/corrs/np_ew.root");  //NP+EW corrections

    int tag = Year.Contains("15") ? 15 : 16;

    TH1D *hEW = (TH1D*) fNPEW->Get(Form("ew%d_ak4_y%d", tag, y));
    TH1D *hNP = (TH1D*) fNPEW->Get(Form("np%d_ak4_y%d", tag, y));

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



void plotRatio()
{
    TFile *fTh    = TFile::Open("theorFiles/cmsJetsNLO.root");  //NLO predictions
    TFile *fD  = TFile::Open(Form("xFitterTables/data%s.root", year.Data()));

    for(int y = 0; y < 5; ++y) {
        TH1D *hStat = (TH1D*) fD->Get(Form("hStat_y%d",y));
        TH1D *hSysUp = (TH1D*) fD->Get(Form("hSysUp_y%d",y));
        TH1D *hSysDn = (TH1D*) fD->Get(Form("hSysDn_y%d",y));
        TGraphAsymmErrors *gSys = getBand(hStat, hSysUp, hSysDn);

        TString tag = (year == "15") ? "histOld" : "histNew";

        TH1D *hTh = (TH1D*) fTh->Get(tag+Form("CT14Scl_Cnt_y%d",y));
        TH1D *hThScU = (TH1D*) fTh->Get(tag+Form("CT14Scl_Up_y%d",y));
        TH1D *hThScD = (TH1D*) fTh->Get(tag+Form("CT14Scl_Dn_y%d",y));
        TH1D *hThPdfU = (TH1D*) fTh->Get(tag+Form("CT14PDF_Up_y%d",y));
        TH1D *hThPdfD = (TH1D*) fTh->Get(tag+Form("CT14PDF_Dn_y%d",y));

        applyNPEW(hTh,    y, year);
        applyNPEW(hThScU, y, year);
        applyNPEW(hThScD, y, year);
        applyNPEW(hThPdfU, y, year);
        applyNPEW(hThPdfD, y, year);





        hStat->Divide(hTh); //normalize to theory
        hThScU->Divide(hTh); //normalize to theory
        hThScD->Divide(hTh); //normalize to theory
        hThPdfU->Divide(hTh); //normalize to theory
        hThPdfD->Divide(hTh); //normalize to theory

        TCanvas *can = new TCanvas(rn(), "", 600, 400);
        SetTopBottom(0.1, 0.15);
        gStyle->SetOptStat(0);
        can->SetLogx();
        can->SetTicky(1);

        hStat->Draw("axis");

        gSys->SetFillColor(kOrange);
        gSys->Draw("le2 same");


        hThScU->SetLineColor(kRed);
        hThScD->SetLineColor(kRed);
        hThPdfU->SetLineColor(kRed);
        hThPdfD->SetLineColor(kRed);
        hThPdfU->SetLineStyle(2);
        hThPdfD->SetLineStyle(2);



        hThScU->Draw("hist same ][");
        hThScD->Draw("hist same ][");
        hThPdfU->Draw("hist same ][");
        hThPdfD->Draw("hist same ][");

        hStat->SetMarkerStyle(20);
        hStat->SetLineColor(kBlack);
        hStat->SetMarkerColor(kBlack);
        hStat->Draw("e0  same");

        GetXaxis()->SetTitle("Jet p_{T} (GeV)");
        GetXaxis()->SetNoExponent();
        GetXaxis()->SetMoreLogLabels();
        GetYaxis()->SetTitle("Ratio to NLOJet++ CT14");
        GetYaxis()->SetRangeUser(0.1, 2.5);

        SetFTO({20}, {10}, {1.1, 2.1, 0.5, 3.1});


        can->Print(Form("plots/data%s_y%d.pdf", year.Data(), y));
    }



}

void plotJets()
{
   plotRatio(); 



}
