#include <vector>

R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)

//#include "tools.h"

#include "plottingHelper.h"
#include "RemoveOverlaps.h"
#include <fstream>

using namespace PlottingHelper;//pollute the namespace!

TString rn() {return Form("%d",rand());}

vector<TString> yBins = {"|y| < 0.5",  "0.5 < |y| < 1.0",  "1.0 < |y| < 1.5", "1.5 < |y| < 2.0", "2.0 < |y| < 2.5"};

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






TObject *getSafe (TFile *f, TString str)
{
    TObject *obj = f->Get(str);
    if(!obj) {
        cout << "Histo " << str << " not loaded from file " << f->GetName() << endl;
    }
    return obj->Clone(rn());
    //return obj;
}

pair<TGraphAsymmErrors*, TH1D*> getEnvelope (vector<TH1D*> hs)
{
    TH1D *hCnt = (TH1D*) hs[0]->Clone(rn());
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
        gr->SetPoint(i-1, x, y);
        gr->SetPointError(i-1, xerr, xerr, yerr, yerr);
        hCnt->SetBinContent(i, y);
        hCnt->SetBinError(i, 1e-8);
    }
    return make_pair(gr, hCnt);
}


void printEnvelope(vector<vector<TH1D*>> hs, TString fName)
{
    vector<int> bins = {56  , 74  , 97  , 133 , 174 , 220 , 272 , 330 , 395 , 468 , 548 , 638 , 737 , 846 , 967 , 1101, 1248, 1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103, 3450, 3832,};
    vector<int> edges= {3832, 3450, 3103, 2238, 2000};

    ofstream outFile(fName);

    for(int iy = 0; iy < hs[0].size(); ++iy) {
        for(int ipt = 0; ipt < bins.size()-1; ++ipt) {
            int pt1 = bins[ipt];
            int pt2 = bins[ipt+1];
            if(pt2 > edges[iy]) continue;

            double Max = -1e10, Min = 1e10;
            for(int i = 1; i <= hs[0][0]->GetNbinsX(); ++i) {

                double ptL   = hs[0][0]->GetXaxis()->GetBinLowEdge(i);
                double ptH   = hs[0][0]->GetXaxis()->GetBinUpEdge(i);
                if(int(round(ptL)) != pt1) continue;
                if(int(round(ptH)) != pt2) continue;
                for(int n = 0; n < hs.size(); ++n) {
                    auto &h = hs[n][iy];
                    Max = max(Max, h->GetBinContent(i) + h->GetBinError(i));
                    Min = min(Min, h->GetBinContent(i) - h->GetBinError(i));
                }
            }

            double corr    = (Max + Min)/2;
            double corrErr = (Max - Min)/2;

            if(Max == -1e10 || Min == 1e10) {corr = 1; corrErr = 0; }
            //cout << iy*0.5 <<" "<< (iy+1)*0.5 <<" "<< pt1 <<" "<< pt2 << " "<<  corr <<" "<< corrErr << " : " << Min <<" "<< Max << endl;
            double relUnc = (corr != 0) ? corrErr / corr * 100 : 0;
            outFile << setw(10) << iy*0.5 <<" "<< setw(10) << (iy+1)*0.5 <<" "<<setw(10)<< pt1 <<" "<<setw(10)<< pt2 << " "<<setw(12)<<  corr <<" "<<setw(12)<< relUnc << " "<<setw(12)<< -relUnc << endl;
        }

    }
    outFile.close();
}





//Load correction from file fHad / fNoHad
pair<vector<TH1D*>, vector<TH1D*>> loadCorr(TFile *fHad, TFile *fNoHad)
{
    //contains nameIndex + y bin
    vector<TH1D*> NPak4, NPak7;
    for(int R : vector<int>({4, 7})) {
        //cout << R << endl;
        vector<TH1D*> NPnow;
        for(int y = 0; y < 5; ++y) {
            TH1D *hHad   = (TH1D*) getSafe(fHad,Form("CMS_2019_incJets/ak%d_y%d",R,y));
            TH1D *hNoHad = (TH1D*) getSafe(fNoHad,Form("CMS_2019_incJets/ak%d_y%d",R,y));
            if(!hNoHad) {
                cout << "Radek " << hNoHad << endl;
                exit(1);
            }

            /*
            //Rebin to the smaller histogram
            if(hHad->GetNbinsX() > hNoHad->GetNbinsX())
                hHad = rebinCorr(hHad, hNoHad);
            if(hHad->GetNbinsX() < hNoHad->GetNbinsX())
                hNoHad = rebinCorr(hNoHad, hHad);
            */

            //cout <<"Radek " <<  hHad->GetNbinsX() << " "<< hNoHad->GetNbinsX() << endl;
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


    void loadData ()
    {
        vector<TString> names = {"QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8",  "QCD_Pt-15to7000_TuneCP1_Flat_13TeV_pythia8", "QCD_Pt-15to7000_Tune4C_Flat_13TeV_pythia8",
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

    void plotNP (int R)
    {

        TCanvas *can = new TCanvas(rn(), "", 1200, 500);
        SetTopBottom(0.05, 0.14);
        SetLeftRight(0.10, 0.05);

        auto &hNP    = (R == 4) ? NPak4 : NPak7;
        auto &hNPmpi = (R == 4) ? NPak4Hg_MPI : NPak7Hg_MPI;

        DividePad({1,1,1,1,1}, {1});

        for(int y = 0; y < 5; ++y) {
            can->cd(y+1);
            gPad->SetLogx();
            gPad->SetTicks(1,1);
            gStyle->SetOptStat(0);
            hNP[0][y]->SetLineColor(kBlack);
            hNP[0][y]->SetLineWidth(2);
            hNP[0][y]->Draw("axis");


            TH1D *hCnt;
            TGraphAsymmErrors *gr;
            tie(gr,hCnt)= getEnvelope({hNP[0][y],hNP[1][y],hNP[2][y],hNP[3][y]} );
            //TGraphAsymmErrors *gr = getEnvelope({hNP[0][y],hNP[1][y], hNP[2][y], hNP[3][y], hNP[4][y]}  );
            gr->SetFillStyle(1001);
            gr->SetFillColor(kOrange);
            gr->SetLineColor(kOrange+2);
            gr->Draw("2 same");
            hCnt->SetLineColor(kOrange+2);
            hCnt->SetMarkerColor(kOrange+2);
            hCnt->Draw("same");


            vector<int> Cols = {kBlue, kMagenta, kGreen+2, kCyan+2};
            for(int i = 1; i < hNP.size(); ++i) {
                hNP[i][y]->SetLineColor(Cols[i-1]);
                hNP[i][y]->SetMarkerColor(Cols[i-1]);
                hNP[i][y]->SetLineWidth(2);
                hNP[i][y]->Draw("hist e ][ same");
                hNP[i][y]->Draw("hist  ][ same");
            }


            hNPmpi[y]->SetLineColor(kRed);
            hNPmpi[y]->SetMarkerColor(kRed);
            hNPmpi[y]->SetLineWidth(2);
            hNPmpi[y]->SetLineStyle(2);
            hNPmpi[y]->Draw("hist e ][ same");
            hNPmpi[y]->Draw("hist  ][ same");

            /*
            hNP[3][y]->SetLineColor(kRed);
            hNP[3][y]->Draw("hist e  ][ same");
            hNP[3][y]->Draw("hist   ][ same");
            */

            hNP[0][y]->Draw("hist e ][ same");
            hNP[0][y]->Draw("hist  ][ same");


            GetYaxis()->SetRangeUser(0.83, 1.34);
            GetXaxis()->SetMoreLogLabels();
            GetXaxis()->SetNoExponent();


            if(y == 0) GetXaxis()->SetRangeUser(74, 3400);
            if(y == 1) GetXaxis()->SetRangeUser(74, 3103);
            if(y == 2) GetXaxis()->SetRangeUser(74, 2941);
            if(y == 3) GetXaxis()->SetRangeUser(74, 2366);
            if(y == 4) GetXaxis()->SetRangeUser(74, 1410);


            if(y == 0) GetYaxis()->SetTitle("Non-perturbative correction factor");
            if(y == 4) GetXaxis()->SetTitle("Jet p_{T} (GeV)");


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
                leg->AddEntry(hNP[1][y], "Py8 CP1", "l");
                leg->AddEntry(hNP[2][y], "Py8 4C", "l");
                DrawLegends({leg}, true);
            }


            else if(y == 2) {
                auto leg = newLegend(kPos7);
                leg->SetTextSize(PxFontToRel(20));
                leg->AddEntry((TObject*)0, "", "h");
                leg->AddEntry(hNP[3][y], "Hg++ CUETHS1", "l");
                leg->AddEntry(hNP[4][y], "Hg7 CH3", "l");
                leg->AddEntry(hNPmpi[y], "Hg7 CH3 (Had)", "l");
                DrawLegends({leg}, true);
            }
            else if(y == 3) {
                auto leg = newLegend(kPos7);
                leg->SetTextSize(PxFontToRel(20));
                leg->AddEntry((TObject*)0, "", "h");
                leg->AddEntry(gr, "unc. envelope", "fl");
                leg->AddEntry((TObject*)0, "+central val.", "");
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

void plotNP ()
{

    plotter pl;
    pl.loadData();

    //printEnvelope( {pl.NPak7[0],pl.NPak7[1],pl.NPak7[2],pl.NPak7[3]}, "np16_ak7.txt");
    //printEnvelope( {pl.NPak4[0],pl.NPak4[1],pl.NPak4[2],pl.NPak4[3]}, "np16_ak4.txt");
    //printEnvelope(pl.NPak4, "np16_ak4.txt");

    pl.plotNP(7);
    pl.plotNP(4);




}
