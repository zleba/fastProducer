//********************************************************************
//     
//     fnlo-tk-h1diffpdf.cc
//     Program to read fastNLO v2 tables and derive
//     QCD cross sections using PDFs e.g. from LHAPDF
//     
//********************************************************************
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cfloat>
//#include "fastnlotk/fastNLODiffReader.h"
//#include "fastNLODiffAlphas.h"
#include "fastnlotk/fastNLOLHAPDF.h"

#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"

#include "plottingHelper.h"
using namespace PlottingHelper;


// Function prototype for flexible-scale function 
double Function_Mu(double s1, double s2 );

double Function_q(double q, double pt );
double Function_pt(double q, double pt );


TString rn() {return Form("%d",rand());}

using namespace std;

vector<TH1D*> readHisto(fastNLOLHAPDF &fnlo)
{
    fnlo.SetScaleFactorsMuRMuF(1.0, 1.0);
    fnlo.CalcCrossSection();
    vector<double> xs = fnlo.GetCrossSection();
    
    map<double, vector<double>> bins, xSec;
    for(int k = 0; k < xs.size(); ++k) {
        //cout << "Helenka " << k <<" "<< fnlodiff.GetObsBinLoBound(k,0) << endl;
        double etaDn = fnlo.GetObsBinLoBound(k,0);
        double ptDn  = fnlo.GetObsBinLoBound(k,1);
        double etaUp = fnlo.GetObsBinUpBound(k,0);
        double ptUp  = fnlo.GetObsBinUpBound(k,1);

        bins[etaDn].push_back(ptDn);
        if(k == xs.size() - 1 || fnlo.GetObsBinLoBound(k,0) != fnlo.GetObsBinLoBound(k+1,0))
            bins[etaDn].push_back(ptUp);

        xSec[etaDn].push_back(xs[k]);
        //cout << etaAvg <<" "<<ptAvg << " "<< xs[0][k] <<" "<< xs[1][k] <<" "<< xs[2][k] <<" "<< xs[2][k]<<  endl;
    }
    cout << "First part done" << endl;

    vector<TH1D*> hists;
    for(auto obj : bins) {
        double etaDn = obj.first;
        vector<double> binning = obj.second;

        TH1D * h = new TH1D(rn(), Form("%g", etaDn), binning.size()-1, binning.data());

        cout << etaDn << endl;

        for(int i = 0; i < xSec.at(etaDn).size(); ++i) {
            h->SetBinContent(i+1, xSec.at(etaDn)[i]);
            h->SetBinError(i+1, 0);
        }

        hists.push_back(h);
    }
    cout << "Second part done" << endl;

    return hists;

}

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

void printHisto(TH1D *h)
{
    for(int i = 1; i < h->GetNbinsX(); ++i) {
        cout << i << " "<<  h->GetBinLowEdge(i) <<" "<<  h->GetXaxis()->GetBinUpEdge(i) <<" "<< h->GetBinContent(i) <<  endl;
    }
}


//__________________________________________________________________________________________________________________________________


vector<TString> yBins = {"|y| < 0.5",  "0.5 < |y| < 1",  "1 < |y| < 1.5", "1.5 < |y| < 2", "2 < |y| < 2.5"};


int main(int argc, char** argv){
	

  // namespaces
  using namespace std;
  using namespace say;		// namespace for 'speaker.h'-verbosity levels
  using namespace fastNLO;	// namespace for fastNLO constants

	SetGlobalVerbosity(DEBUG);


    //vector<TH1D*> readHisto(fastNLOLHAPDF &fnlo, TString tabName, TString pdfName)

    //say::SetGlobalVerbosity(say::DEBUG);
    fastNLOLHAPDF fnloOld("theorFiles/fastnlo-cms-incjets-arxiv-1605.04436-xsec001.tab", "CT14nlo", 0);
    fastNLOLHAPDF fnloNew("theorFiles/InclusiveNJets_fnl5362h_v23_fix.tab", "CT14nlo", 0);
    fnloOld.SetContributionON(fastNLO::kFixedOrder,0,true);
    fnloOld.SetContributionON(fastNLO::kFixedOrder,1,true);
    fnloOld.SetUnits(fastNLO::kPublicationUnits);

    fnloNew.SetContributionON(fastNLO::kFixedOrder,0,true);
    fnloNew.SetContributionON(fastNLO::kFixedOrder,1,true);
    fnloNew.SetUnits(fastNLO::kPublicationUnits);


    int nMem = fnloOld.GetNPDFMembers();


    vector<TH1D*> oldH = readHisto(fnloOld);
    vector<TH1D*> newH = readHisto(fnloNew);


    TCanvas *can = new TCanvas(rn(), "", 1000, 350);
    gStyle->SetOptStat(0);
    SetLeftRight(0.07, 0.03);
    SetTopBottom(0.05, 0.13);

    DividePad( {1,1,1,1,1}, {1});

    for(int yB = 0; yB < 5; ++yB) {

        TH1D *oldHH = rebin(oldH[yB], newH[yB]);


        oldHH->Divide(newH[yB]);

        can->cd(yB+1);
        gPad->SetLogx();
        oldHH->Draw("][");
        GetYaxis()->SetRangeUser(0.99, 1.01);

        if(yB == 4)
            GetXaxis()->SetTitle("Jet p_{T} (GeV)");
        GetYaxis()->SetTitle("#sigma^{NLO2015} / #sigma^{NLO2016} ");
        GetXaxis()->SetRangeUser(100, 2600);
        GetXaxis()->SetNoExponent();
        GetXaxis()->SetMoreLogLabels();
        SetFTO({16}, {10}, {1.25, 2.3, 0.3, 4.2});


        TLegend *leg = newLegend(kPos8c);
        leg->AddEntry((TObject*)nullptr, yBins[yB], "h");
        if(yB == 0) {
            leg->AddEntry((TObject*)nullptr, "AK4 jets", "h");
            leg->AddEntry(oldHH, "NLOJET++", "l");
            leg->AddEntry((TObject*)nullptr, "from fastNLO", "");
        }
        DrawLegends({leg}, false);

        TLine *line = new TLine;
        line->SetLineStyle(2);
        line->DrawLine(100, 1.001, 2600, 1.001);
        line->DrawLine(100, 0.999, 2600, 0.999);


        for(int i = 1; i < oldHH->GetNbinsX(); ++i) {
            cout << i << " "<<  oldHH->GetBinLowEdge(i) <<" "<<  oldHH->GetBinContent(i) <<" "<< newH[yB]->GetBinContent(i) <<  endl;
        }

    }
    can->SaveAs("tableCheck.pdf");



  
  return 0;
}



//__________________________________________________________________________________________________________________________________

double Function_q(double q, double pt ){
  // --- fastNLO user: This is an example function
  //     to demonstrate how you might perform the
  //     definition of the scales using a 
    return sqrt(q*q);
}
double Function_pt(double q, double pt ){
  // --- fastNLO user: This is an example function
  //     to demonstrate how you might perform the
  //     definition of the scales using a 
    return 1.1*sqrt(pt*pt);
}



double Function_Mu(double s1, double s2 ){
  // --- fastNLO user: This is an example function
  //     to demonstrate how you might perform the
  //     definition of the scales using a 
  //     'flexible-scale'-table
   //double mu = s1*exp(0.3*s2);

   //double mu = sqrt(s1*s1/4. + s2*s2);
   double mu = sqrt(s1*s1 + s2*s2);
   return mu;
}

//__________________________________________________________________________________________________________________________________
