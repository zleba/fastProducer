R__ADD_INCLUDE_PATH($PlH_DIR/PlottingHelper)
R__LOAD_LIBRARY($PlH_DIR/plottingHelper_C.so)

#include "tools.h"

#include "plottingHelper.h"
using namespace PlottingHelper;//pollute the namespace!



TString rn() {return Form("%d",rand());}



struct Fitter {
    map<TString, vector<vector<TGraph*>> > grY, grPt;
    map<TString, vector<TGraph*> > grAll;
    void load(TString fName) {
        TFile *f = TFile::Open(fName);
        //vector<TString> pdfNames = {"CT14nlo", "HERAPDF20_NLO", "NNPDF31_nnlo"};
        vector<TString> pdfNames = {"CT14nnlo"};//, "NNPDF31_nnlo"};
        for(auto nPdf : pdfNames) {
            vector<vector<TGraph*>> gyVec(7), gptVec(7);
            vector<TGraph*> gVecAll(7);
            for(int s = 0; s < 7; ++s) {
                gyVec[s].resize(4);

                for(int y = 0; y < gyVec[s].size(); ++y) {
                    gyVec[s][y] = (TGraph*) f->Get(nPdf+Form("_Y%d_scale%d_all", y, s));
                    if(!gyVec[s][y]) {
                        cout << "Histogram cannot be read " << nPdf+Form("_Y%d_scale%d", y, s) << endl;
                        exit(1);
                    }
                }

                gptVec[s].resize(ptBinsAs.size()-1);

                /*
                for(int ipt = 0; ipt < gptVec[s].size(); ++ipt) {
                    gptVec[s][ipt] = (TGraph*) f->Get(nPdf+Form("_Pt%d_scale%d", ipt, s));
                    if(!gptVec[s][ipt]) {
                        cout << "Histogram cannot be read " << nPdf+Form("_Pt%d_scale%d", ipt, s) << endl;
                        exit(1);
                    }
                }
                */


                gVecAll[s] = (TGraph*) f->Get(nPdf+Form("_scale%d", s));
            }
            grY[nPdf] = gyVec;
            //grPt[nPdf] = gptVec;
            grAll[nPdf] = gVecAll;
        }
    }

    vector<double> fitAs(TGraph *gr, TF1 **Fit = nullptr)
    {
        auto &fit = *Fit;
        fit = new TF1(rn(), "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x", -15, 10);
        gr->Fit(fit, "N", "", 0.110, 0.220);
        double asMin = fit->GetMinimumX();
        double chi2Min = fit->Eval(asMin);
        double shLow  = fit->GetX(chi2Min + 1, -15, asMin);
        double shHigh = fit->GetX(chi2Min + 1, asMin, 10);

        double errL = asMin - shLow;
        double errH = shHigh - asMin;
        //cout << "Helenka min " << asMin << " "<< errL <<" "<< errH << endl;
        return {asMin, errH, errL};
    }

    vector<double> fitAsScale(TString pdfName, TString type,  int var)
    {
        vector<TGraph*> grNow;
        if(type == "tot") {
            for(int s = 0; s < 7; ++s)
                grNow.push_back(grAll.at(pdfName)[s]);
        }
        else if(type == "y") {
            for(int s = 0; s < 7; ++s)
                grNow.push_back(grY.at(pdfName)[s][var]);
        }
        else if(type == "pt") {
            for(int s = 0; s < 7; ++s)
                grNow.push_back(grPt.at(pdfName)[s][var]);
        }
        else {
            cout << "Something wrong - Radek" << endl;
            exit(1);
        }


        TF1 *fD;
        vector<double> res = fitAs(grNow[0], &fD);
        double scaleMin = 100, scaleMax = -100;
        for(int s = 0; s < 7; ++s) {
            TF1 *fDd;
            vector<double> resS = fitAs(grNow[s], &fDd);
            scaleMin = min(scaleMin, resS[0]);
            scaleMax = max(scaleMax, resS[0]);
        }
        return {res[0], res[1], res[2],  scaleMax - res[0], res[0] - scaleMin };
    }


    void plotChi2Fit(TString pdfName, int y) {

       cout << __LINE__ << endl;

        TCanvas *can = new TCanvas(rn(), "", 600, 600);
        can->SetLeftMargin(0.15);


       cout << __LINE__ << endl;
        vector<TGraph*> grNow;
        for(int s = 0; s < 7; ++s)
            grNow.push_back(grY.at(pdfName)[s][y]);

       cout << __LINE__ << endl;
        grNow[0]->SetMarkerStyle(20);

       cout << __LINE__ << endl;
        grNow[0]->Draw("ap");
        for(int s = 1; s < 7; ++s) {

            TF1 *fNow = (TF1*) grNow[s]->GetListOfFunctions()->At(0);
            fNow->SetLineColor(kBlack);
            fNow->SetLineWidth(1);
            fNow->SetLineStyle(2);
            grNow[s]->SetMarkerStyle(4);
            grNow[s]->Draw("p same");

        }

       cout << __LINE__ << endl;

        TF1 *fRes = nullptr;
        vector<double> res = fitAs(grNow[0], &fRes);
        cout << "Pusinka " << fRes << endl;
        double val = fRes->Eval(res[0]);
        TLine *line = new TLine;
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->DrawLine(res[0]-res[2], val+1, res[0] + res[1], val+1);
        line->DrawLine(res[0], val, res[0], val+1);


       cout << __LINE__ << endl;

        double scaleMin = res[0], scaleMax = res[0];
        for(int s = 1; s < 7; ++s) {
            TF1 *f;
            vector<double> resS = fitAs(grNow[s], &f);
            scaleMin = min(scaleMin, resS[0]);
            scaleMax = max(scaleMax, resS[0]);
        }
        line->SetLineColor(kBlue);
        cout << "ScaleMin ScaleMax " << scaleMin <<" "<< scaleMax << endl;
        line->DrawLine(scaleMin, val+0, scaleMax, val+0);


       cout << __LINE__ << endl;

        GetXaxis()->SetTitle("#alpha_{S}(M_{Z})");
        GetYaxis()->SetTitle("#chi^{2}");

       cout << __LINE__ << endl;
        TLegend *leg = newLegend(kPos8c);
        leg->AddEntry((TObject*)0, "|y| < 0.5", "h");
        leg->AddEntry((TObject*)0, Form("#alpha_{s} = %.4f^{+%.4f}_{-%.4f} (exp)^{+%.4f}_{-%.4f}(scl)", res[0], res[1], res[2], scaleMax-res[0], res[0]-scaleMin), "h");
        leg->SetFillStyle(1001);
        DrawLegends({leg});

        can->SaveAs("asFitChi2_" + pdfName + ".pdf");
    }


    void plotYdep(TString pdfName)
    {
        TGraphAsymmErrors *gr = new TGraphAsymmErrors();
        TGraphAsymmErrors *grS = new TGraphAsymmErrors();
        for(int y = 0; y < 4; ++y) {
            //vector<double> res = fitAs(grY.at(pdfName)[0][y]);

            vector<double> res = fitAsScale(pdfName, "y", y);

            cout << y <<" "<<res[0] << endl;
            gr->SetPoint(y, y*0.5 + 0.25, res[0]);
            gr->SetPointError(y, 0.25, 0.25, res[2], res[1]);

            grS->SetPoint(y, y*0.5 + 0.25, res[0]);
            grS->SetPointError(y, 0.25, 0.25, hypot(res[2],res[4]) , hypot(res[1],res[3]) );
            //grS->SetPointError(y, 0.25, 0.25, res[4] , res[3] );

        }

        vector<double> res = fitAsScale(pdfName, "tot",  -1);
        cout << "Radek " << res[0] << endl;
        TGraphAsymmErrors *grTot = new TGraphAsymmErrors();
        grTot->SetPoint(0, 1.0, res[0]);
        grTot->SetPointError(0, 1.0, 1.0, res[2], res[1]);
        TGraphAsymmErrors *grTotS = new TGraphAsymmErrors();
        grTotS->SetPoint(0, 1.0, res[0]);
        //grTotS->SetPointError(0, 1.0, 1.0, res[4], res[3]);
        grTotS->SetPointError(0, 1.0, 1.0, hypot(res[2],res[4]) , hypot(res[1],res[3]) );


        TCanvas *can = new TCanvas(rn(), "", 600, 600);
        can->SetLeftMargin(0.15);


        //grTotS->SetFillColor(kBlue);
        grTotS->SetLineColor(kRed);
        grTotS->SetLineWidth(0);
        //grTotS->SetFillStyle(3315);
        //grTotS->Draw("a le2");
        grTotS->Draw("a le2");

        grTot->SetFillColor(kOrange);
        grTot->SetFillStyle(3144);
        grTot->SetLineColor(kRed);
        grTot->SetLineWidth(2);
        //grTot->Draw("le2 same");
        /*
        TLine *l = new TLine();
        l->SetLineColor(kRed);
        l->DrawLine(0, res[0], 2, res[0]);
        */


        //grS->SetLineColor(kRed);
        grS->Draw("p same");

        gr->SetMarkerStyle(20);
        gr->Draw("p same");

        GetXaxis()->SetRangeUser(0, 2);
        GetYaxis()->SetRangeUser(0.105, 0.125);

        GetXaxis()->SetTitle("|y|");
        GetYaxis()->SetTitle("#alpha^{NLO}_{S}(M_{Z})");

        UpdateFrame();

        can->SaveAs("asYdep_" + pdfName + ".pdf");

    }


    void plotPtDep(TString pdfName)
    {
        TGraphAsymmErrors *gr = new TGraphAsymmErrors();
        TGraphAsymmErrors *grS = new TGraphAsymmErrors();
        for(int ipt = 0; ipt < ptBinsAs.size()-1; ++ipt) {
            //vector<double> res = fitAs(grY.at(pdfName)[0][y]);

            //cout << "Before " << ipt << endl;
            vector<double> res = fitAsScale(pdfName, "pt", ipt);
            //cout << "After " << endl;

            cout << "Radek " << ipt <<" "<<res[0] << endl;

            double ptMin = ptBinsAs[ipt];
            double ptMax = ptBinsAs[ipt+1];
            double ptAvg = sqrt(ptMin*ptMax);

            gr->SetPoint(ipt, ptAvg, res[0]);
            gr->SetPointError(ipt, ptAvg-ptMin, ptMax-ptAvg, res[2], res[1]);

            grS->SetPoint(ipt, ptAvg, res[0]);
            grS->SetPointError(ipt, ptAvg-ptMin, ptMax-ptAvg, res[4], res[3]);
        }

        vector<double> res = fitAsScale(pdfName, "tot",  -1);
        cout << "Radek " << res[0] << endl;
        TGraphAsymmErrors *grTot = new TGraphAsymmErrors();

        double ptMin = ptBinsAs[0];
        double ptMax = ptBinsAs[ptBinsAs.size()-1];
        double ptAvg = sqrt(ptMin*ptMax);


        grTot->SetPoint(0, ptAvg, res[0]);
        grTot->SetPointError(0, ptAvg-ptMin, ptMax-ptAvg, res[2], res[1]);
        TGraphAsymmErrors *grTotS = new TGraphAsymmErrors();
        grTotS->SetPoint(0, ptAvg, res[0]);
        grTotS->SetPointError(0, ptAvg-ptMin, ptMax-ptAvg, res[4], res[3]);

        grTot->SetFillColor(kOrange);
        grTot->SetLineColor(kRed);
        grTot->Draw("a le2");

        grTotS->SetFillColor(kBlue);
        grTotS->SetLineColor(kRed);
        grTotS->Draw("le2 same");


        gr->Draw("* same");
        grS->SetLineColor(kRed);
        grS->Draw("* same");

        GetXaxis()->SetRangeUser(97, 3103);
        GetYaxis()->SetRangeUser(0.095, 0.125);

        GetXaxis()->SetTitle("|y|");
        GetYaxis()->SetTitle("#alpha^{NLO}_{S}(M_{Z})");

        UpdateFrame();


    }



};

void plotChi2()
{
    Fitter fitter;
    fitter.load("chi2.root");
    //fitter.load("chi2noCorr.root");
    fitter.plotYdep("CT14nnlo");
    //fitter.plotYdep("NNPDF31_nnlo");

    return;
    fitter.plotChi2Fit("CT14nnlo", 0);


    //fitter.plotPtDep("CT14nlo");

    //fitter.plotYdep("HERAPDF20_NLO");



}
