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

    TH1D *hEW = dynamic_cast<TH1D*>( fNPEW->Get(Form("ew%d_ak4_y%d", tag, y)));
    TH1D *hNP = dynamic_cast<TH1D*>( fNPEW->Get(Form("np%d_ak4_y%d", tag, y)));
    if(!hEW || !hNP) {
        cout << "Histogram is missing in np_ew.root file" << endl;
        exit(0);
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
        cout << "Histogram not found :" << Tag + Form("_ak4_y%d",  y) << endl;
        exit(0);
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




vector<TString> yBins = {"|y| < 0.5",  "0.5 < |y| < 1",  "1 < |y| < 1.5", "1.5 < |y| < 2", "2 < |y| < 2.5"};

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

        TH1D *hThNNLO = (TH1D*) hTh->Clone(rn());
        TH1D *hThNLL  = (TH1D*) hTh->Clone(rn());

        applyKfactor(hThNLL, y, "kFactorNLL");
        applyKfactor(hThNNLO, y, "kFactorNNLO");



        hStat->Divide(hTh); //normalize to theory
        hThScU->Divide(hTh); //normalize to theory
        hThScD->Divide(hTh); //normalize to theory
        hThPdfU->Divide(hTh); //normalize to theory
        hThPdfD->Divide(hTh); //normalize to theory

        hThNLL->Divide(hTh);
        hThNNLO->Divide(hTh);


        TCanvas *can = new TCanvas(rn(), "", 550, 400);
        SetTopBottom(0.1, 0.15);
        gStyle->SetOptStat(0);
        can->SetLogx();
        can->SetTicky(1);

        hStat->Draw("axis");

        gSys->SetFillColor(kOrange);
        gSys->SetLineColor(kBlack);
        gSys->SetLineStyle(9);
        gSys->SetLineWidth(2);
        gSys->Draw("le2 same");


        hThScU->SetLineColor(kRed);
        hThScD->SetLineColor(kRed);
        hThPdfU->SetLineColor(kRed);
        hThPdfD->SetLineColor(kRed);
        hThPdfU->SetLineStyle(2);
        hThPdfD->SetLineStyle(2);

        hThNLL->SetLineWidth(2);
        hThNNLO->SetLineWidth(2);
        hThNLL->SetLineColor(kBlue);
        hThNNLO->SetLineColor(kMagenta);


        hThScU->Draw("hist same ][");
        hThScD->Draw("hist same ][");
        hThPdfU->Draw("hist same ][");
        hThPdfD->Draw("hist same ][");

        hThNLL->Draw("hist same ][");
        hThNNLO->Draw("hist same ][");



        hStat->SetMarkerStyle(20);
        hStat->SetLineColor(kBlack);
        hStat->SetMarkerColor(kBlack);
        hStat->Draw("e0  same");

        GetXaxis()->SetTitle("Jet p_{T} (GeV)");
        GetXaxis()->SetNoExponent();
        GetXaxis()->SetMoreLogLabels();
        GetYaxis()->SetTitle("Ratio to NLOJet++ CT14");
        GetYaxis()->SetRangeUser(0.1, 2.5);

        SetFTO({20}, {10}, {1.15, 2.1, 0.3, 2.73});



        auto leg = newLegend(kPos7);

        leg->SetNColumns(2);
        //leg->SetMargin (0.4);

        leg->AddEntry((TObject*)nullptr, "Inclusive jets R = 0.4", "");    leg->AddEntry(hThScU,  "NLO scale unc.", "l");
        leg->AddEntry((TObject*)nullptr, yBins[y], "");                    leg->AddEntry(hThPdfU, "NLO PDF unc.", "l");

        leg->AddEntry(hStat, "Data + (stat unc.)", "pe");                  leg->AddEntry(hThNLL,  "NLO+NLL", "l");
        leg->AddEntry(gSys , "Exp. unc", "f");                              leg->AddEntry(hThNNLO, "NNLO", "l");

        DrawLegends({leg}, true);

        UpdateFrame();

        can->Print(Form("plots/data%s_y%d.pdf", year.Data(), y));
    }
}


void compareRatio(TString year0, TString year1)
{
    TFile *fTh = TFile::Open("theorFiles/cmsJetsNLO.root");  //NLO predictions
    TFile *fD  = TFile::Open(Form("xFitterTables/data%s.root", year0.Data()));
    TFile *fD1 = TFile::Open(Form("xFitterTables/data%s.root", year1.Data()));

    for(int y = 0; y < 5; ++y) {
        TH1D *hStat = (TH1D*) fD->Get(Form("hStat_y%d",y));
        TH1D *hStat1= (TH1D*) fD1->Get(Form("hStat_y%d",y));
        TH1D *hSysUp = (TH1D*) fD->Get(Form("hSysUp_y%d",y));
        TH1D *hSysDn = (TH1D*) fD->Get(Form("hSysDn_y%d",y));
        TGraphAsymmErrors *gSys = getBand(hStat, hSysUp, hSysDn);

        TString tag0 = (year0 == "15") ? "histOld" : "histNew";

        TH1D *hTh = (TH1D*) fTh->Get(tag0+Form("CT14Scl_Cnt_y%d",y));
        TH1D *hThScU = (TH1D*) fTh->Get(tag0+Form("CT14Scl_Up_y%d",y));
        TH1D *hThScD = (TH1D*) fTh->Get(tag0+Form("CT14Scl_Dn_y%d",y));
        TH1D *hThPdfU = (TH1D*) fTh->Get(tag0+Form("CT14PDF_Up_y%d",y));
        TH1D *hThPdfD = (TH1D*) fTh->Get(tag0+Form("CT14PDF_Dn_y%d",y));

        TString tag1 = (year1 == "15") ? "histOld" : "histNew";
        TH1D *hTh1 = (TH1D*) fTh->Get(tag1+Form("CT14Scl_Cnt_y%d",y))->Clone(rn());

        applyNPEW(hTh,    y, year0);
        applyNPEW(hThScU, y, year0);
        applyNPEW(hThScD, y, year0);
        applyNPEW(hThPdfU, y, year0);
        applyNPEW(hThPdfD, y, year0);
        applyNPEW(hTh1,    y, year1);

        TH1D *hThNNLO = (TH1D*) hTh->Clone(rn());
        TH1D *hThNLL  = (TH1D*) hTh->Clone(rn());

        applyKfactor(hThNLL, y, "kFactorNLL");
        applyKfactor(hThNNLO, y, "kFactorNNLO");



        hStat->Divide(hTh); //normalize to theory
        hStat1->Divide(hTh1); //normalize to theory
        hThScU->Divide(hTh); //normalize to theory
        hThScD->Divide(hTh); //normalize to theory
        hThPdfU->Divide(hTh); //normalize to theory
        hThPdfD->Divide(hTh); //normalize to theory

        hThNLL->Divide(hTh);
        hThNNLO->Divide(hTh);


        TCanvas *can = new TCanvas(rn(), "", 550, 400);
        SetTopBottom(0.1, 0.15);
        gStyle->SetOptStat(0);
        can->SetLogx();
        can->SetTicky(1);

        hStat->Draw("axis");

        gSys->SetFillColor(kOrange);
        gSys->SetLineColor(kBlack);
        gSys->SetLineStyle(9);
        gSys->SetLineWidth(2);
        gSys->Draw("le2 same");


        hThScU->SetLineColor(kRed);
        hThScD->SetLineColor(kRed);
        hThPdfU->SetLineColor(kRed);
        hThPdfD->SetLineColor(kRed);
        hThPdfU->SetLineStyle(2);
        hThPdfD->SetLineStyle(2);

        hThNLL->SetLineWidth(2);
        hThNNLO->SetLineWidth(2);
        hThNLL->SetLineColor(kBlue);
        hThNNLO->SetLineColor(kMagenta);


        hThScU->Draw("hist same ][");
        hThScD->Draw("hist same ][");
        hThPdfU->Draw("hist same ][");
        hThPdfD->Draw("hist same ][");

        hThNLL->Draw("hist same ][");
        hThNNLO->Draw("hist same ][");



        hStat->SetMarkerStyle(20);
        hStat->SetLineColor(kBlack);
        hStat->SetMarkerColor(kBlack);
        hStat->Draw("e0  same");

        hStat1->SetMarkerStyle(20);
        hStat1->SetLineColor(kBlue);
        hStat1->SetMarkerColor(kBlue);
        hStat1->Draw("e0  same");


        GetXaxis()->SetTitle("Jet p_{T} (GeV)");
        GetXaxis()->SetNoExponent();
        GetXaxis()->SetMoreLogLabels();
        GetYaxis()->SetTitle("Ratio to NLOJet++ CT14");
        GetYaxis()->SetRangeUser(0.1, 2.5);

        if(year0 == "16" && year1 == "16e" || true) {
            if(y == 0) GetXaxis()->SetRangeUser(97, 3103);
            if(y == 1) GetXaxis()->SetRangeUser(97, 2940);
            if(y == 2) GetXaxis()->SetRangeUser(97, 2787);
            if(y == 3) GetXaxis()->SetRangeUser(97, 2000);
        }

        SetFTO({20}, {10}, {1.15, 2.1, 0.3, 2.73});



        auto leg = newLegend(kPos7);

        leg->SetNColumns(2);
        //leg->SetMargin (0.4);

        leg->AddEntry((TObject*)nullptr, "Inclusive jets R = 0.4", "");    leg->AddEntry(hThScU,  "NLO scale unc.", "l");
        leg->AddEntry((TObject*)nullptr, yBins[y], "");                    leg->AddEntry(hThPdfU, "NLO PDF unc.", "l");

        leg->AddEntry(hStat,  "Data "+year0+" + (stat unc.)", "pe");       leg->AddEntry(hThNLL,  "NLO+NLL", "l");
        leg->AddEntry(hStat1, "Data "+year1+" + (stat unc.)", "pe");       leg->AddEntry(hThNNLO, "NNLO", "l");

        leg->AddEntry(gSys , "Exp. unc", "f");

        DrawLegends({leg}, true);

        UpdateFrame();

        can->Print(Form("plots/data%s_y%d.pdf", (year0+year1).Data(), y));
    }
}




void plotJets()
{
   //plotRatio(); 
   //compareRatio("15", "16");
   
   //compareRatio("16", "15");
   compareRatio("16", "16e");

}
