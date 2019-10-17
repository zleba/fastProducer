R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)

#include "tools.h"

#include "plottingHelper.h"
#include "RemoveOverlaps.h"
using namespace PlottingHelper;//pollute the namespace!

TString rn() {return Form("%d",rand());}

TString year = "16";


vector<TString> yBins = {"|y| < 0.5",  "0.5 < |y| < 1",  "1 < |y| < 1.5", "1.5 < |y| < 2", "2 < |y| < 2.5"};

void myAssertFun(bool cond,  int line)
{
    if(!cond) {
        cout << "Issue at line " << line << endl;
        exit(1);
    }
}

#define myAssert(cond) myAssertFun(cond,   __LINE__)



TH1D *rebin(TH1D *h, TH1D *hTemp) //second is template
{
    TH1D *hNew = (TH1D *) hTemp->Clone(rn());
    hNew->Reset();

    for(int i = 1; i <= h->GetNbinsX(); ++i) {
        double xCnt = h->GetBinCenter(i);

        double iBin = hNew->FindBin(xCnt);

        if(iBin == 0 || iBin == hNew->GetNbinsX()+1)
            continue;

        double v = h->GetBinContent(i);
        double w = h->GetBinWidth(i);

        double orgVal = hNew->GetBinContent(iBin);
        hNew->SetBinContent(iBin, orgVal + v*w);
    }

    hNew->Scale(1, "width");
    return hNew;
}



//Modify histogram h
void removeEmpty(TH1D *h, TH1D *hTemp)
{
    for(int i = 1; i <= h->GetNbinsX(); ++i) { 
        if(hTemp->GetBinContent(i) < 1e-14)
            h->SetBinContent(i, 0);
    }
}


TGraphAsymmErrors *getBandTot(TH1D *hNom, vector<TH1D*> hSh)
{
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(hNom->GetNbinsX());
    for(int i = 1; i <= hNom->GetNbinsX(); ++i) {
        double cnt = hNom->GetBinContent(i);

        double dn = 0;
        double up = 0;
        for(auto h : hSh) {
            double v = h->GetBinContent(i);
            up =  hypot(up, max(0.,v - cnt));
            dn =  hypot(dn, max(0.,cnt - v));
        }
        gr->SetPoint(i-1, hNom->GetBinCenter(i), cnt);
        gr->SetPointError(i-1, hNom->GetBinWidth(i)/2, hNom->GetBinWidth(i)/2, dn, up);
    }
    return gr;
}


void plotTheorUnc(TString Tag, TString pdf)
{
    TFile *fTh;
    
    if(Tag == "16ak7")
        fTh = TFile::Open("theorFiles/cmsJetsNLO_AK7.root");  //NLO predictions
    else if(Tag == "16ak4")
        fTh = TFile::Open("theorFiles/cmsJetsNLO_AK4.root");  //NLO predictions
    else {
        cout << "Bad file" << endl;
        exit(1);
    }

    TFile *fD  = TFile::Open(Form("xFitterTables/data%s.root", Tag.Data()));

    //int yMax = (Tag != "16ak7") ? 5 : 4;
    int yMax =  4;

    //vector<TH1D*> hStat(yMax);

    TCanvas *can = new TCanvas(rn(), "", 650, 300);
    SetTopBottom(0.1, 0.15);
    DividePad({1,1,1,1}, {1});
    gStyle->SetOptStat(0);


    for(int y = 0; y < yMax; ++y) {
        //TString tag = (Year == "15") ? "histOld" : "histNew";

        //can->cd(y+1);

        TH1D *hStat = (TH1D*) fD->Get(Form("hStat_y%d",y));
        myAssert(hStat);

        TH1D *hTh = (TH1D*) fTh->Get("hist"+pdf+Form("_Scl_Cnt_y%d",y));
        myAssert(hTh);
        TH1D *hThScU  = (TH1D*) fTh->Get("hist"+pdf+Form("_Scl_Up_y%d",y));
        TH1D *hThScD  = (TH1D*) fTh->Get("hist"+pdf+Form("_Scl_Dn_y%d",y));
        TH1D *hThPdfU = (TH1D*) fTh->Get("hist"+pdf+Form("_PDF_Up_y%d",y));
        TH1D *hThPdfD = (TH1D*) fTh->Get("hist"+pdf+Form("_PDF_Dn_y%d",y));
        TH1D *hThAsU    = (TH1D*) fTh->Get("hist"+pdf+Form("_As_Up_y%d",y));
        myAssert(hThAsU);
        TH1D *hThAsD    = (TH1D*) fTh->Get("hist"+pdf+Form("_As_Dn_y%d",y));
        myAssert(hThAsD);

        hTh = rebin(hTh, hStat);
        hThScU = rebin(hThScU, hStat);
        hThScD = rebin(hThScD, hStat);
        hThPdfU = rebin(hThPdfU, hStat);
        hThPdfD = rebin(hThPdfD, hStat);
        hThAsU = rebin(hThAsU, hStat);
        hThAsD = rebin(hThAsD, hStat);

        //hAsU->Print("all");
        //hAsD->Print("all");
        //exit(0);



        //hStat->Divide(hTh); //normalize to theory
        hThScU->Divide(hTh); //normalize to theory
        hThScD->Divide(hTh); //normalize to theory
        hThPdfU->Divide(hTh); //normalize to theory
        hThPdfD->Divide(hTh); //normalize to theory
        hThAsU->Divide(hTh); //normalize to theory
        hThAsD->Divide(hTh); //normalize to theory

        hTh->Divide(hTh);


        //hThNLL->Divide(hTh);
        //hThNNLO->Divide(hTh);


        can->cd(y+1);

        gPad->SetLogx();
        gPad->SetTicky(1);

        hStat->Draw("axis");

        TGraphAsymmErrors *grTot = getBandTot(hTh, {hThScU, hThScD, hThPdfU, hThPdfD, hThAsU, hThAsD});
        grTot->SetFillColor(kOrange);
        grTot->Draw("e2 same");


        //gSys->SetFillColor(kOrange);
        //gSys->SetLineColor(kBlack);
        //gSys->SetLineStyle(9);
        //gSys->SetLineWidth(2);
        //gSys->Draw("le2 same");


        hThScU->SetLineColor(kRed);
        hThScD->SetLineColor(kRed);
        hThPdfU->SetLineColor(kRed);
        hThPdfD->SetLineColor(kRed);
        hThPdfU->SetLineStyle(2);
        hThPdfD->SetLineStyle(2);

        hThAsU->SetLineColor(kBlue);
        hThAsD->SetLineColor(kBlue);
        hThAsU->SetLineStyle(1);
        hThAsD->SetLineStyle(1);

        //hThNLL->SetLineWidth(2);
        //hThNNLO->SetLineWidth(2);
        //hThNLL->SetLineColor(kBlue);
        //hThNNLO->SetLineColor(kMagenta);


        hThScU->Draw("hist same ][");
        hThScD->Draw("hist same ][");
        hThPdfU->Draw("hist same ][");
        hThPdfD->Draw("hist same ][");
        hThAsU->Draw("hist same ][");
        hThAsD->Draw("hist same ][");




        if(y == 3) GetXaxis()->SetTitle("Jet p_{T} (GeV)");
        GetXaxis()->SetNoExponent();
        GetXaxis()->SetMoreLogLabels();
        GetYaxis()->SetTitle("Theoretical Uncertainty");
        //GetYaxis()->SetRangeUser(0.1, 2.5);
        GetYaxis()->SetRangeUser(0.7, 1.5);

        SetFTO({16}, {10}, {1.15, 2.1, 0.3, 2.73});

        if(y == 0) GetXaxis()->SetRangeUser(97, 3103);
        if(y == 1) GetXaxis()->SetRangeUser(97, 2940);
        if(y == 2) GetXaxis()->SetRangeUser(97, 2787);
        if(y == 3) GetXaxis()->SetRangeUser(97, 2000);


        if(y == 0) {
            auto leg = newLegend(kPos7);
            leg->AddEntry((TObject*)0,  "", "h");
            int R = Tag.Contains("ak4") ? 4 : 7;
            leg->AddEntry((TObject*)0, "#sqrt{s} = 13TeV" , "h");
            leg->AddEntry((TObject*)0, Form("anti-k_{T} (R=0.%d)", R) , "h");
            leg->AddEntry((TObject*)0,  pdf+" PDF", "h");
            DrawLegends({leg}, true);
        }



        if(y == 1) {
            auto leg = newLegend(kPos7);
            leg->AddEntry((TObject*)0,  "", "h");
            leg->AddEntry(hThScU,  "Scale unc.", "l");
            leg->AddEntry(hThPdfU, "PDF unc.", "l");
            leg->AddEntry(hThAsU,  "#alpha_{S} unc.", "l");
            leg->AddEntry(grTot,  "Total unc.", "f");
            DrawLegends({leg}, true);
        }


        DrawLatexUp(-1,  yBins[y]);

        RemoveOverlaps(gPad, GetXaxis());

        /*
        auto leg = newLegend(kPos7);

        leg->SetNColumns(2);
        //leg->SetMargin (0.4);

        double R = (Year != "16ak7") ? 0.4 : 0.7;

        leg->AddEntry((TObject*)nullptr, Form("Inclusive jets R = %g", R), "");    leg->AddEntry(hThScU,  "NLO scale unc.", "l");
        leg->AddEntry((TObject*)nullptr, yBins[y], "");                       leg->AddEntry(hThPdfU, "NLO PDF unc.", "l");

        leg->AddEntry(hStat, "Data + (stat unc.)", "pe");   
        if( R != 0.7) leg->AddEntry(hThNLL,  "NLO+NLL", "l");
        else leg->AddEntry((TObject*)0,  "", "");
        leg->AddEntry(gSys , "Exp. unc", "f"); 
        if( R != 0.7) leg->AddEntry(hThNNLO, "NNLO", "l");
        else leg->AddEntry((TObject*)0,  "", "");

        DrawLegends({leg}, true);
        */

        UpdateFrame();

    }
    can->Print(Form("plots/Theor_%s_%s.pdf", pdf.Data(), Tag.Data()));
}




void plotJetsLog(TString Year)
{
    TFile *fTh;
    
    if(Year == "16ak7")
        fTh = TFile::Open("theorFiles/cmsJetsNLO_AK7.root");  //NLO predictions
    else
        fTh = TFile::Open("theorFiles/cmsJetsNLO.root");  //NLO predictions

    TFile *fD  = TFile::Open(Form("xFitterTables/data%s.root", Year.Data()));

    int yMax = (Year != "16ak7") ? 5 : 4;
    yMax = 4;


    TCanvas *can = new TCanvas(rn(), "", 600, 600);
    can->SetLogx();
    can->SetLogy();
    can->SetTicky(1);
    SetTopBottom(0.07, 0.13);
    SetLeftRight(0.15, 0.1);
    gStyle->SetOptStat(0);


    //TH1D *hSysUp = (TH1D*) fD->Get(Form("hSysUp_y%d",y));
    //TH1D *hSysDn = (TH1D*) fD->Get(Form("hSysDn_y%d",y));
    //TGraphAsymmErrors *gSys = getBand(hStat, hSysUp, hSysDn);


    //vector<double> yShifts = {1, 1e2, 1e4, 1e6, 1e8};
    vector<double> yShifts = {1, 1e-1, 1e-2, 1e-3, 1e-4};

    vector<TH1D*> hStat(yMax);

    TH1D *hThLeg;
    for(int y = 0; y < yMax; ++y) {
        hStat[y] = (TH1D*) fD->Get(Form("hStat_y%d",y));

        TString tag = (Year == "15") ? "histOld" : "histNew";

        TH1D *hTh = (TH1D*) fTh->Get(tag+Form("CT14Scl_Cnt_y%d",y));
        TH1D *hThScU = (TH1D*) fTh->Get(tag+Form("CT14Scl_Up_y%d",y));
        TH1D *hThScD = (TH1D*) fTh->Get(tag+Form("CT14Scl_Dn_y%d",y));
        TH1D *hThPdfU = (TH1D*) fTh->Get(tag+Form("CT14PDF_Up_y%d",y));
        TH1D *hThPdfD = (TH1D*) fTh->Get(tag+Form("CT14PDF_Dn_y%d",y));

        hTh = rebin(hTh, hStat[y]);
        hThScU = rebin(hThScU, hStat[y]);
        hThScD = rebin(hThScD, hStat[y]);
        hThPdfU = rebin(hThPdfU, hStat[y]);
        hThPdfD = rebin(hThPdfD, hStat[y]);


        //Move by factors


        applyNPEW(hTh,    y, Year);
        applyNPEW(hThScU, y, Year);
        applyNPEW(hThScD, y, Year);
        applyNPEW(hThPdfU, y, Year);
        applyNPEW(hThPdfD, y, Year);

        TH1D *hThNNLO = (TH1D*) hTh->Clone(rn());
        TH1D *hThNLL  = (TH1D*) hTh->Clone(rn());

        if(Year != "16ak7") {
            applyKfactor(hThNLL, y,  "kFactorNLL");
            applyKfactor(hThNNLO, y, "kFactorNNLO");
        }

        /*
        hStat->Divide(hTh); //normalize to theory
        hThScU->Divide(hTh); //normalize to theory
        hThScD->Divide(hTh); //normalize to theory
        hThPdfU->Divide(hTh); //normalize to theory
        hThPdfD->Divide(hTh); //normalize to theory

        hThNLL->Divide(hTh);
        hThNNLO->Divide(hTh);
        */

        hStat[y]->Scale(yShifts[y]);
        hTh->Scale(yShifts[y]);


        if(y == 0) hStat[y]->Draw("axis");

        hTh->SetLineColor(kRed);

        removeEmpty(hTh, hStat[y]);

        hTh->Draw("same ][");

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

        /*
        hThScU->Draw("hist same ][");
        hThScD->Draw("hist same ][");
        hThPdfU->Draw("hist same ][");
        hThPdfD->Draw("hist same ][");


        if(Year != "16ak7") {
            hThNLL->Draw("hist same ][");
            hThNNLO->Draw("hist same ][");
        }
        */

        hStat[y]->SetMarkerStyle(20+y);
        hStat[y]->SetLineColor(kBlack);
        hStat[y]->SetMarkerColor(kBlack);
        hStat[y]->Draw("e0  same");

        hThLeg = hTh;
    }

    GetXaxis()->SetTitle("Jet  p_{T} (GeV)");
    GetXaxis()->SetNoExponent();
    GetXaxis()->SetMoreLogLabels();
    GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}dy (pb/GeV)");
    //GetYaxis()->SetRangeUser(0.1, 2.5);
    GetYaxis()->SetRangeUser(1e-9, 1e6);

    SetFTO({20}, {10}, {1.15, 2.1, 0.3, 2.83});


    GetXaxis()->SetRangeUser(97, 3403);

    DrawLatexUp(0.8, "#bf{<35.9 fb^{-1} (13 TeV)}", 27, "r");
    DrawLatexUp(-1, "    #bf{CMS}", 27, "l");

    UpdateFrame();


    auto leg = newLegend(kPos9);

    //leg->SetNColumns(2);
    //leg->SetMargin (0.4);

    double R = (Year != "16ak7") ? 0.4 : 0.7;

    leg->AddEntry((TObject*)nullptr, Form("Anti-k_{t} jets (R = %g)", R), "");
    //leg->AddEntry(hThScU,  "NLO scale unc.", "l");
    //leg->AddEntry((TObject*)nullptr, yBins[y], "");
    //leg->AddEntry(hThPdfU, "NLO PDF unc.", "l");

    vector<TString> yNames = {"(x10^{0})", "(x10^{-1})","(x10^{-2})", "(x10^{-3})", "(x10^{-4})"};
    for(int y = 0; y < yMax; ++y)
        leg->AddEntry(hStat[y], yBins[y]+" "+yNames[y], "pe");   

    leg->AddEntry(hThLeg, "NLO CT14", "l");   
    /*
    if( R != 0.7) leg->AddEntry(hThNLL,  "NLO+NLL", "l");
    else leg->AddEntry((TObject*)0,  "", "");
    leg->AddEntry(gSys , "Exp. unc", "f"); 
    if( R != 0.7) leg->AddEntry(hThNNLO, "NNLO", "l");
    else leg->AddEntry((TObject*)0,  "", "");
    */

    DrawLegends({leg}, true);


    GetXaxis()->SetTitleSize(PxFontToRel(27));
    GetYaxis()->SetTitleSize(PxFontToRel(27));



    can->Print(Form("plots/dataLog%s.pdf", Year.Data()));

}









void plotRatio(TString Year)
{
    TFile *fTh;
    
    if(Year == "16ak7")
        fTh = TFile::Open("theorFiles/cmsJetsNLO_AK7.root");  //NLO predictions
    else
        fTh = TFile::Open("theorFiles/cmsJetsNLO.root");  //NLO predictions

    TFile *fD  = TFile::Open(Form("xFitterTables/data%s.root", Year.Data()));

    int yMax = (Year != "16ak7") ? 5 : 4;

    for(int y = 0; y < yMax; ++y) {
        TH1D *hStat = (TH1D*) fD->Get(Form("hStat_y%d",y));
        TH1D *hSysUp = (TH1D*) fD->Get(Form("hSysUp_y%d",y));
        TH1D *hSysDn = (TH1D*) fD->Get(Form("hSysDn_y%d",y));
        TGraphAsymmErrors *gSys = getBand(hStat, hSysUp, hSysDn);

        TString tag = (Year == "15") ? "histOld" : "histNew";

        TH1D *hTh = (TH1D*) fTh->Get(tag+Form("CT14Scl_Cnt_y%d",y));
        TH1D *hThScU = (TH1D*) fTh->Get(tag+Form("CT14Scl_Up_y%d",y));
        TH1D *hThScD = (TH1D*) fTh->Get(tag+Form("CT14Scl_Dn_y%d",y));
        TH1D *hThPdfU = (TH1D*) fTh->Get(tag+Form("CT14PDF_Up_y%d",y));
        TH1D *hThPdfD = (TH1D*) fTh->Get(tag+Form("CT14PDF_Dn_y%d",y));

        hTh = rebin(hTh, hStat);
        hThScU = rebin(hThScU, hStat);
        hThScD = rebin(hThScD, hStat);
        hThPdfU = rebin(hThPdfU, hStat);
        hThPdfD = rebin(hThPdfD, hStat);


        applyNPEW(hTh,    y, Year);
        applyNPEW(hThScU, y, Year);
        applyNPEW(hThScD, y, Year);
        applyNPEW(hThPdfU, y, Year);
        applyNPEW(hThPdfD, y, Year);

        TH1D *hThNNLO = (TH1D*) hTh->Clone(rn());
        TH1D *hThNLL  = (TH1D*) hTh->Clone(rn());

        if(Year != "16ak7") {
            applyKfactor(hThNLL, y,  "kFactorNLL");
            applyKfactor(hThNNLO, y, "kFactorNNLO");
        }


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

        if(Year != "16ak7") {
            hThNLL->Draw("hist same ][");
            hThNNLO->Draw("hist same ][");
        }



        hStat->SetMarkerStyle(20);
        hStat->SetLineColor(kBlack);
        hStat->SetMarkerColor(kBlack);
        hStat->Draw("e0  same");

        GetXaxis()->SetTitle("Jet p_{T} (GeV)");
        GetXaxis()->SetNoExponent();
        GetXaxis()->SetMoreLogLabels();
        GetYaxis()->SetTitle("Ratio to NLOJet++ CT14");
        //GetYaxis()->SetRangeUser(0.1, 2.5);
        GetYaxis()->SetRangeUser(0.7, 1.5);

        SetFTO({20}, {10}, {1.15, 2.1, 0.3, 2.73});

        //const double magNum = 97; 
        const double magNum = 74; 

        if(y == 0) GetXaxis()->SetRangeUser(magNum, 3103);
        if(y == 1) GetXaxis()->SetRangeUser(magNum, 2940);
        if(y == 2) GetXaxis()->SetRangeUser(magNum, 2787);
        if(y == 3) GetXaxis()->SetRangeUser(magNum, 2000);






        auto leg = newLegend(kPos7);

        leg->SetNColumns(2);
        //leg->SetMargin (0.4);

        double R = (Year != "16ak7") ? 0.4 : 0.7;

        leg->AddEntry((TObject*)nullptr, Form("Inclusive jets R = %g", R), "");    leg->AddEntry(hThScU,  "NLO scale unc.", "l");
        leg->AddEntry((TObject*)nullptr, yBins[y], "");                       leg->AddEntry(hThPdfU, "NLO PDF unc.", "l");

        leg->AddEntry(hStat, "Data + (stat unc.)", "pe");   
        if( R != 0.7) leg->AddEntry(hThNLL,  "NLO+NLL", "l");
        else leg->AddEntry((TObject*)0,  "", "");
        leg->AddEntry(gSys , "Exp. unc", "f"); 
        if( R != 0.7) leg->AddEntry(hThNNLO, "NNLO", "l");
        else leg->AddEntry((TObject*)0,  "", "");

        DrawLegends({leg}, true);

        UpdateFrame();

        can->Print(Form("plots/data%s_y%d.pdf", Year.Data(), y));
    }
}


void compareRatio(TString year0, TString year1)
{
    TFile *fTh = nullptr;
    
    if(year0.Contains("ak7"))
        fTh = TFile::Open("theorFiles/cmsJetsNLO_AK7.root");
    else
        fTh = TFile::Open("theorFiles/cmsJetsNLO.root");  //NLO predictions



    TFile *fD  = TFile::Open(Form("xFitterTables/data%s.root", year0.Data()));
    TFile *fD1 = TFile::Open(Form("xFitterTables/data%s.root", year1.Data()));

    cout << "Fname " << Form("xFitterTables/data%s.root", year0.Data()) << endl;

    for(int y = 0; y < 4; ++y) {
        TH1D *hStat = (TH1D*) fD->Get(Form("hStat_y%d",y));
        //hStat->Print("all");
        TH1D *hStat1= (TH1D*) fD1->Get(Form("hStat_y%d",y));
        TH1D *hSysUp = (TH1D*) fD->Get(Form("hSysUp_y%d",y));
        TH1D *hSysDn = (TH1D*) fD->Get(Form("hSysDn_y%d",y));
        TGraphAsymmErrors *gSys = getBand(hStat, hSysUp, hSysDn);

        TString tag0 = year0.Contains("15") ? "histOld" : "histNew";

        TH1D *hTh = (TH1D*) fTh->Get(tag0+Form("CT14Scl_Cnt_y%d",y));
        TH1D *hThScU = (TH1D*) fTh->Get(tag0+Form("CT14Scl_Up_y%d",y));
        TH1D *hThScD = (TH1D*) fTh->Get(tag0+Form("CT14Scl_Dn_y%d",y));
        TH1D *hThPdfU = (TH1D*) fTh->Get(tag0+Form("CT14PDF_Up_y%d",y));
        TH1D *hThPdfD = (TH1D*) fTh->Get(tag0+Form("CT14PDF_Dn_y%d",y));


        hTh = rebin(hTh, hStat);
        hThScU = rebin(hThScU, hStat);
        hThScD = rebin(hThScD, hStat);
        hThPdfU = rebin(hThPdfU, hStat);
        hThPdfD = rebin(hThPdfD, hStat);


        TString tag1 =  year1.Contains("15")  ? "histOld" : "histNew";
        cout << "Reading " << y << endl;
        TH1D *hTh1 = (TH1D*) fTh->Get(tag1+Form("CT14Scl_Cnt_y%d",y))->Clone(rn());

        applyNPEW(hTh,    y, year0);
        applyNPEW(hThScU, y, year0);
        applyNPEW(hThScD, y, year0);
        applyNPEW(hThPdfU, y, year0);
        applyNPEW(hThPdfD, y, year0);
        applyNPEW(hTh1,    y, year1);

        TH1D *hThNNLO = (TH1D*) hTh->Clone(rn());
        TH1D *hThNLL  = (TH1D*) hTh->Clone(rn());

        if(!year0.Contains("ak7")) {
            applyKfactor(hThNLL, y, "kFactorNLL");
            applyKfactor(hThNNLO, y, "kFactorNNLO");
        }

        //hStat->Print("all");
        //hTh->Print("all");

        //hTh1 = rebin(hTh1, hStat1);


        hStat->Divide(hTh); //normalize to theory
        hStat1->Divide(hTh1); //normalize to theory
        hThScU->Divide(hTh); //normalize to theory
        hThScD->Divide(hTh); //normalize to theory
        hThPdfU->Divide(hTh); //normalize to theory
        hThPdfD->Divide(hTh); //normalize to theory


        hStat->Print("all");
        hStat1->Print("all");
        //hTh1->Print("all");
        //hStat1->Print("all");

        //exit(0);



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

        hStat1->SetMarkerStyle(kOpenStar);
        hStat1->SetLineColor(kBlack);
        hStat1->SetMarkerColor(kBlack);
        //hStat1->Scale(0.8);
        hStat1->Draw("e0  same");

        //hStat1->Print("all"); //Helenka printing

        GetXaxis()->SetTitle("Jet p_{T} (GeV)");
        GetXaxis()->SetNoExponent();
        GetXaxis()->SetMoreLogLabels();
        GetYaxis()->SetTitle("Ratio to NLOJet++ CT14");
        GetYaxis()->SetRangeUser(0.5, 2.5);

        if((year0 == "16" && year1 == "16e") || true) {
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

void plotAsScan(TString pdfName)
{
    TFile *fTh = TFile::Open("cmsJetsAsScan.root");  //NLO predictions
    TFile *fD  = TFile::Open(Form("xFitterTables/data%s.root", year.Data()));

    for(int y = 0; y < 5; ++y) {
        TH1D *hStat = (TH1D*) fD->Get(Form("hStat_y%d",y));
        TH1D *hSysUp = (TH1D*) fD->Get(Form("hSysUp_y%d",y));
        TH1D *hSysDn = (TH1D*) fD->Get(Form("hSysDn_y%d",y));
        TGraphAsymmErrors *gSys = getBand(hStat, hSysUp, hSysDn);

	cout << "In the middle " <<y << endl;
        vector<TH1D*> vTh;
        //for(int as = 111; as <= 123; ++as) {
        int idCnt = 0;
        for(const double & as : pdfAsVals.at(pdfName)) {
            int asI = round(as*1000);
            if(asI < 118) ++idCnt;
            vTh.push_back( (TH1D*) fTh->Get(pdfName+ Form("_y%d_as0%d_scale0_pdf0",y,asI) ));
        }
	cout << "In the middle2 " <<y << endl;
        //cout << "Id cnt is " << idCnt << endl;
        //exit(0);

        for( auto &h : vTh) {
            applyNPEW(h,  y, year);
            applyKfactor(h, y, "kFactorNLL");
            //applyKfactor(h, y, "kFactorNNLO");
        }
	cout << "In the middle3 " <<y << endl;

        //applyNPEW(hTh,    y, year);
        //applyNPEW(hThScU, y, year);
        //applyNPEW(hThScD, y, year);
        //applyNPEW(hThPdfU, y, year);
        //applyNPEW(hThPdfD, y, year);

        //applyKfactor(hThNLL, y, "kFactorNLL");
        //applyKfactor(hThNNLO, y, "kFactorNNLO");



        hStat->Divide(vTh[idCnt]); //normalize to theory

        /*
        if(y == 0) {
            hStat->Print("all");
            //vTh[7]->Print("all");
            exit(0);
        }
        */


        for(int i = 0; i < vTh.size(); ++i) {
            if(i == idCnt) continue;
            vTh[i]->Divide(vTh[idCnt]);
        }

        //hThNLL->Divide(hTh);
        //hThNNLO->Divide(hTh);


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



        for(int i = 0; i < vTh.size(); ++i) {
            if(i == idCnt) continue;
            vTh[i]->Draw("hist same ][");
        }


        hStat->SetMarkerStyle(20);
        hStat->SetLineColor(kBlack);
        hStat->SetMarkerColor(kBlack);
        hStat->Draw("e0  same");

        GetXaxis()->SetTitle("Jet p_{T} (GeV)");
        GetXaxis()->SetNoExponent();
        GetXaxis()->SetMoreLogLabels();
        GetYaxis()->SetTitle("Ratio to NLOJet++ " + pdfName);
        GetYaxis()->SetRangeUser(0.5, 1.5);

        SetFTO({20}, {10}, {1.15, 2.1, 0.3, 2.73});


        if(y == 0) GetXaxis()->SetRangeUser(97, 3103);
        else if(y == 1) GetXaxis()->SetRangeUser(97, 2940);
        else if(y == 2) GetXaxis()->SetRangeUser(97, 2787);
        else if(y == 3) GetXaxis()->SetRangeUser(97, 2000);



        double asMin = pdfAsVals.at(pdfName)[0];
        double asMax = pdfAsVals.at(pdfName)[pdfAsVals.at(pdfName).size()-1];

        auto leg = newLegend(kPos7);

        leg->SetNColumns(2);
        //leg->SetMargin (0.4);

        leg->AddEntry((TObject*)nullptr, "Inclusive jets R = 0.4", "");    leg->AddEntry(vTh[idCnt],  "NLO+NLL+NP+EW", "l");
        leg->AddEntry((TObject*)nullptr, yBins[y], "");                    leg->AddEntry(vTh[idCnt], Form("#alpha_{S} = %g - %g", asMin, asMax), "l");

        leg->AddEntry(hStat, "Data + (stat unc.)", "pe");                  leg->AddEntry((TObject*)nullptr,  "", "");
        leg->AddEntry(gSys , "Exp. unc", "f");                             leg->AddEntry((TObject*)nullptr, "", "");

        DrawLegends({leg}, true);


        UpdateFrame();

        can->Print(Form("plots/dataScan%s_y%d_%s.pdf", year.Data(), y, pdfName.Data()));
    }
}


//Plot theor unc for several PDFs
void plotTheorUncAll()
{
    for(auto pdf :  {"CT14nlo", "HERAPDF20_NLO", "NNPDF31_nlo", "ABMP16_5_nlo"}) {
        plotTheorUnc("16ak4", pdf);
        plotTheorUnc("16ak7", pdf);
    }

}



void plotJets()
{
   //plotRatio("16"); 
   //compareRatio("16", "15"); 
   //compareRatio("16", "16m");
   //compareRatio("16ak7", "15ak7");
   //plotRatio("16"); 

    plotTheorUncAll();
    //plotTheorUnc("16ak4", "CT14nlo");
    //plotTheorUnc("16ak7", "CT14nlo");


   //plotJetsLog("16n");
   //plotJetsLog("16");

   return;

     plotAsScan("HERAPDF20_NLO");
     plotAsScan("NNPDF31_nnlo");
     plotAsScan("CT14nlo");
   //compareRatio("16", "16e");

}
