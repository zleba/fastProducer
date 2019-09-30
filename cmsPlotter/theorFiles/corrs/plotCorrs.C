R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)

#include "plottingHelper.h"
using namespace PlottingHelper;//pollute the namespace!

TString rn() {return Form("%d",rand());}

TString year = "16";


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


struct corrPlotter {

    vector<TH1D*> ew15, ew16, np15, np16;

    void init() {
        TFile *file = TFile::Open("np_ew.root");

        ew15.resize(5, nullptr);
        ew16.resize(5, nullptr);
        np15.resize(5, nullptr);
        np16.resize(5, nullptr);


        for(int y = 0; y < 5; ++y) {
            ew15[y] = dynamic_cast<TH1D*>(file->Get(Form("ew15_ak4_y%d", y)));
            ew16[y] = dynamic_cast<TH1D*>(file->Get(Form("ew16_ak4_y%d", y)));
            np15[y] = dynamic_cast<TH1D*>(file->Get(Form("np15_ak4_y%d", y)));
            np16[y] = dynamic_cast<TH1D*>(file->Get(Form("np16_ak4_y%d", y)));
        }
    }

    vector<TString> yBins = {"|y| < 0.5",  "0.5 < |y| < 1",  "1 < |y| < 1.5", "1.5 < |y| < 2", "2 < |y| < 2.5"};


    void plotYearComp(TString tag)
    {
        gStyle->SetOptStat(0);
        TCanvas *can = new TCanvas(rn(), "",  1000, 350);
        SetLeftRight(0.05, 0.05);
        SetTopBottom(0.05, 0.13);

        DividePad({1,1,1,1,1}, {1});

        for(int y = 0; y < 5; ++y) {
            TH1D *h1, *h2;
            if(tag == "EW") {
                h1 = ew15[y];
                h2 = ew16[y];
            }
            else {
                h1 = np15[y];
                h2 = np16[y];
            }
            can->cd(y+1);
            gPad->SetLogx();

            h1->SetLineColor(kRed);
            h2->SetLineColor(kBlue);
            h1->SetMarkerColor(kRed);
            h2->SetMarkerColor(kBlue);


            h2->Draw();
            h1->Draw("same ][");

            GetFrame()->SetTitle("");
            GetYaxis()->SetRangeUser(0.93, 1.07);

            if(y == 4)
                GetXaxis()->SetTitle("Jet p_{T} (GeV)");
            GetXaxis()->SetNoExponent();
            GetXaxis()->SetMoreLogLabels();
            GetYaxis()->SetTitle("Correction");

            SetFTO({16}, {10}, {1.25, 2.3, 0.3, 3.3});


            TLegend *leg = newLegend(kPos7);
            leg->AddEntry((TObject*)nullptr, yBins[y], "h");
            if(y == 0) {
                leg->AddEntry((TObject*)nullptr, "AK4 jets", "h");
                leg->AddEntry(h1, tag + " 2015", "l");
                leg->AddEntry(h2, tag + " 2016", "l");
            }
            DrawLegends({leg}, false);
        }
        can->SaveAs("plots/YearByYear"+tag+".pdf");
    }

    void plotEW()
    {
        gStyle->SetOptStat(0);
        TCanvas *can = new TCanvas(rn(), "",  600, 600);
        can->SetTickx(1);
        can->SetTicky(1);

        SetLeftRight(0.15, 0.05);
        SetTopBottom(0.05, 0.13);

        ew16[0]->SetLineColor(kBlack);
        ew16[1]->SetLineColor(kRed);
        ew16[2]->SetLineColor(kBlue);
        ew16[3]->SetLineColor(kMagenta);


        for(int y = 0; y < 4; ++y) {
            can->cd(y+1);
            gPad->SetLogx();


            if(y == 0)
                ew16[y]->Draw("][");
            else
                ew16[y]->Draw("][ same");
        }
        
        GetFrame()->SetTitle("");
        GetYaxis()->SetNdivisions(505);
        GetYaxis()->SetRangeUser(0.93, 1.115);
        GetXaxis()->SetRangeUser(97, 3103);

        GetXaxis()->SetTitle("Jet p_{T} (GeV)");
        GetXaxis()->SetNoExponent();
        GetXaxis()->SetMoreLogLabels();
        GetYaxis()->SetTitle("Electroweak Correction Factor");

        SetFTO({24}, {10}, {1.25, 2.3, 0.3, 3.3});

        DrawLatexUp(0.85, "13 TeV", -1, "r");

        TLegend *leg = newLegend(kPos7);
        //leg->AddEntry((TObject*)nullptr, yBins[y], "h");
        leg->AddEntry((TObject*)nullptr, "Dittmaier, Huss, Speckner", "h");
        leg->AddEntry((TObject*)nullptr, "Anti-k_{T} jets (R=0.4)", "h");

        leg->AddEntry(ew16[0], "|y| < 0.5", "l");
        leg->AddEntry(ew16[1], "0.5 < |y| < 1.0", "l");
        leg->AddEntry(ew16[2], "1.0 < |y| < 1.5", "l");
        leg->AddEntry(ew16[3], "1.5 < |y| < 2.0", "l");


        DrawLegends({leg}, false);

        can->SaveAs("plots/EW_ak4.pdf");
    }







    void plotNPvsEW()
    {
        gStyle->SetOptStat(0);
        TCanvas *can = new TCanvas(rn(), "",  1000, 350);
        SetLeftRight(0.05, 0.05);
        SetTopBottom(0.05, 0.13);

        DividePad({1,1,1,1,1}, {1});

        for(int y = 0; y < 5; ++y) {
            TH1D *ew, *np, *tot;
            ew = ew16[y];
            np = np16[y];
            can->cd(y+1);
            gPad->SetLogx();

            tot = (TH1D*) ew->Clone(rn());
            tot->Multiply(np);


            np->SetLineColor(kRed);
            ew->SetLineColor(kBlue);
            np->SetMarkerColor(kRed);
            ew->SetMarkerColor(kBlue);
            tot->SetLineColor(kBlack);
            tot->SetMarkerColor(kBlack);
            tot->SetLineWidth(2);



            np->Draw("][");
            ew->Draw("same ][");
            tot->Draw("hist same ][");

            GetFrame()->SetTitle("");
            GetYaxis()->SetRangeUser(0.93, 1.07);

            if(y == 4)
                GetXaxis()->SetTitle("Jet p_{T} (GeV)");
            GetXaxis()->SetNoExponent();
            GetXaxis()->SetMoreLogLabels();
            GetYaxis()->SetTitle("Correction");

            SetFTO({16}, {10}, {1.25, 2.3, 0.3, 3.3});




            TLegend *leg = newLegend(kPos7);
            leg->AddEntry((TObject*)nullptr, yBins[y], "h");
            if(y == 0) {
                leg->AddEntry((TObject*)nullptr, "AK4 jets", "h");
                leg->AddEntry(ew, "EW corr.", "l");
                leg->AddEntry(np, "NP corr.", "l");
                leg->AddEntry(tot, "NP+EW corr.", "l");
            }
            DrawLegends({leg}, false);
        }
        can->SaveAs("plots/NPvsEW.pdf");
    }

};




void plotCorrs()
{
    corrPlotter plt;

    plt.init();
    plt.plotEW();

    //plt.plotYearComp("EW");
    //plt.plotYearComp("NP");
    //plt.plotNPvsEW();

}
