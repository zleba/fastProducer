#include <vector>

R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)

//#include "tools.h"

#include "plottingHelper.h"
#include "RemoveOverlaps.h"

using namespace PlottingHelper;//pollute the namespace!

TString rn() {return Form("%d",rand());}

vector<TString> yBins = {"|y| < 0.5",  "0.5 < |y| < 1",  "1 < |y| < 1.5", "1.5 < |y| < 2", "2 < |y| < 2.5"};

using namespace std;


TH1D *rebinCorr(TH1D *h, TH1D *hTemp) //second is template
{
    TH1D *hNew = (TH1D *) hTemp->Clone(rn());
    hNew->Reset();

    for(int i = 1; i <= h->GetNbinsX(); ++i) {
        double xCnt = h->GetBinCenter(i);//from original
        double iBin = hNew->FindBin(xCnt);//find index in new hist

        if(iBin == 0 || iBin == hNew->GetNbinsX()+1)
            continue;

        double v   = h->GetBinContent(i);
        double err = h->GetBinError(i);

        double orgVal = hNew->GetBinContent(iBin);
        if(orgVal != 0) {
            cout << "There is problem "<< endl;
            exit(1);
        }


        hNew->SetBinContent(iBin, v);
        hNew->SetBinError(iBin, err);
    }

    return hNew;
}


















TObject *getSafe(TFile *f, TString str)
{
    TObject *obj = f->Get(str);
    if(!obj) {
        cout << "Histo " << str << " not loaded from file " << f->GetName() << endl;
    }
    return obj->Clone(rn());
    //return obj;
}

TGraphAsymmErrors *getEnvelope(vector<TH1D*> hs)
{
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(hs[0]->GetNbinsX());
    for(int i = 1; i <= hs[0]->GetNbinsX(); ++i) {

        double Max = -1e10, Min = 1e10;
        for(auto h : hs) {
            Max = max(Max, h->GetBinContent(i) + h->GetBinError(i));
            Min = min(Min, h->GetBinContent(i) - h->GetBinError(i));
        }

        double y    = (Max + Min)/2;
        double yerr = (Max - Min)/2;
        double x    = hs[0]->GetBinCenter(i);
        double xerr = hs[0]->GetBinWidth(i)/2;
        gr->SetPoint(i, x, y);
        gr->SetPointError(i, xerr, xerr, yerr, yerr);

    }
    return gr;
}

//Load correction from file fHad / fNoHad
pair<vector<TH1D*>, vector<TH1D*>> loadCorr(TFile *fHad, TFile *fNoHad)
{
    vector<TH1D*> NPak4, NPak7;
    for(int R : vector<int>({4, 7})) {
        //cout << R << endl;
        vector<TH1D*> NPnow;
        for(int y = 0; y < 4; ++y) {
            TH1D *hHad   = (TH1D*) getSafe(fHad,Form("CMS_2019_incJets/ak%d_y%d",R,y));
            TH1D *hNoHad = (TH1D*) getSafe(fNoHad,Form("CMS_2019_incJets/ak%d_y%d",R,y));
            if(!hNoHad) {
                cout << "Radek " << hNoHad << endl;
                exit(1);
            }

            //Rebin to the smaller histogram
            if(hHad->GetNbinsX() > hNoHad->GetNbinsX())
                hHad = rebinCorr(hHad, hNoHad);
            if(hHad->GetNbinsX() < hNoHad->GetNbinsX())
                hNoHad = rebinCorr(hNoHad, hHad);

            cout <<"Radek " <<  hHad->GetNbinsX() << " "<< hNoHad->GetNbinsX() << endl;
            hHad->Divide(hNoHad);
            NPnow.push_back(hHad);
        }
        if(R == 4)  NPak4 = NPnow;
        else        NPak7 = NPnow;
    }
    return make_pair(NPak4, NPak7);
}




struct plotter {
    vector<vector<TH1D*>> NPak4, NPak7;

    vector<TH1D*> NPak4Hg_MPI, NPak7Hg_MPI;


    void loadData() {
        vector<TString> names = {"QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8", // "QCD_Pt-15to7000_TuneCP1_Flat_13TeV_pythia8", "QCD_Pt-15to7000_Tune4C_Flat_13TeV_pythia8",
                                 "QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp",
                                 "QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7" };


        TFile *fMPI = TFile::Open("histos/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7_noHadwMPI.root");
        TFile *fRef = TFile::Open("histos/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7_Had.root");
        tie(NPak4Hg_MPI, NPak7Hg_MPI)  = loadCorr(fRef, fMPI);
        fMPI->Close();
        //fRef->Close();

        //NPak4Hg_MPI[0]->Print("all");
        //exit(0);

        for( auto n : names) {
            TFile *fHad = TFile::Open("histos/"   +n+"_Had.root");
            //cout << "../farm/"+n+"_Had/"+n+"_Had.root" << endl;
            TFile *fNoHad = TFile::Open("histos/"  +n+"_noHad.root");

            auto corr =  loadCorr(fHad, fNoHad);
            NPak4.push_back(corr.first);
            NPak7.push_back(corr.second);
        }
    }

    void plotNP(int R) {

        TCanvas *can = new TCanvas(rn(), "", 1000, 500);
        SetTopBottom(0.05, 0.14);
        SetLeftRight(0.10, 0.05);

        auto &hNP    = (R == 4) ? NPak4 : NPak7;
        auto &hNPmpi = (R == 4) ? NPak4Hg_MPI : NPak7Hg_MPI;

        DividePad({1,1,1,1}, {1});

        for(int y = 0; y < 4; ++y) {
            can->cd(y+1);
            gPad->SetLogx();
            gPad->SetTicks(1,1);
            gStyle->SetOptStat(0);
            hNP[0][y]->SetLineColor(kBlack);
            hNP[0][y]->SetLineWidth(2);
            hNP[0][y]->Draw("axis");


            //TGraphAsymmErrors *gr = getEnvelope({hNP[0][y],hNP[1][y],hNP[2][y],hNP[3][y]} );
            TGraphAsymmErrors *gr = getEnvelope({hNP[0][y],hNP[1][y], hNP[2][y]});
            gr->SetFillStyle(1001);
            gr->SetFillColor(kOrange);
            gr->Draw("2 same");


            hNP[1][y]->Draw("hist e ][ same");
            hNP[1][y]->Draw("hist ][ same");

            hNP[2][y]->SetLineColor(kMagenta);
            hNP[2][y]->Draw("hist e ][ same");
            hNP[2][y]->Draw("hist  ][ same");

            hNPmpi[y]->SetLineColor(kRed);
            hNPmpi[y]->Draw("hist e ][ same");
            hNPmpi[y]->Draw("hist  ][ same");

            /*
            hNP[3][y]->SetLineColor(kRed);
            hNP[3][y]->Draw("hist e  ][ same");
            hNP[3][y]->Draw("hist   ][ same");
            */

            hNP[0][y]->Draw("hist e ][ same");
            hNP[0][y]->Draw("hist  ][ same");


            GetYaxis()->SetRangeUser(0.85, 1.25);
            GetXaxis()->SetMoreLogLabels();
            GetXaxis()->SetNoExponent();

            if(y == 0) GetYaxis()->SetTitle("Non-perturbative correction factor");
            if(y == 3) GetXaxis()->SetTitle("Jet p_{T} (GeV)");


            SetFTO({22}, {9}, {1.5, 2.6, 0.5, 3.7});

            DrawLatexUp(-1.3, yBins[y]);

            if(y == 0) {
                auto leg = newLegend(kPos9);
                leg->AddEntry((TObject*)0, "", "h");
                leg->AddEntry((TObject*)0, "#bf{CMS} #it{simulation}", "h");
                leg->AddEntry((TObject*)0, "#sqrt{s} = 13 TeV", "h");
                leg->AddEntry((TObject*)0, Form("anti-k_{T} (R=0.%d)",R), "h");
                DrawLegends({leg}, true);
            }
            else if(y == 1) {
                auto leg = newLegend(kPos7);
                leg->SetTextSize(PxFontToRel(20));
                leg->AddEntry((TObject*)0, "", "h");
                leg->AddEntry(hNP[0][y], "Py8 CP5", "l");
                //leg->AddEntry(hNP[1][y], "Pythia8 CP1", "l");
                //leg->AddEntry(hNP[2][y], "Pythia8 4C", "l");
                leg->AddEntry(hNP[1][y], "Hg++ CUETHS1", "l");
                leg->AddEntry(hNP[2][y], "Hg7 CH3", "l");
                leg->AddEntry(hNPmpi[y], "Hg7 CH3 (Had)", "l");
                DrawLegends({leg}, true);
            }

            UpdateFrame();

            RemoveOverlaps(gPad, GetXaxis(), true, true);

        }

        can->SaveAs(Form("NP_R%d.pdf", R));

        //NPak7[0][y]->SetLineColor(kRed);
        //NPak7[0][y]->Draw("hist ][ same");


        //gPad->BuildLegend();
    }

};

void plotNP()
{

    plotter pl;
    pl.loadData();
    pl.plotNP(7);
    pl.plotNP(4);






}
