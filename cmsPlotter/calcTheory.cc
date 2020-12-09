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
#include "fastnlotk/fastNLOAlphas.h"

#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"

#include "plottingHelper.h"
#include "tools.h"

using namespace PlottingHelper;


// Function prototype for flexible-scale function 
double Function_Mu(double s1, double s2 );
double Function_q (double q, double pt );
double Function_pt(double q, double pt );


TString rn() {return Form("%d",rand());}

using namespace std;

//All functions to have a list
vector<TH1D*> readHisto(fastNLOAlphas &fnlo);
vector<vector<TH1D*>> getScaleuncHistos(fastNLOAlphas &fnlo);
vector<vector<TH1D*>> getPDFuncHistos(fastNLOAlphas &fnlo);
//TH1D *rebin(TH1D *h, TH1D *hTemp);
//void printHisto(TH1D *h);
void SaveHistos(vector<vector<TH1D*>> hist,  TString tag);
vector<vector<vector<TH1D*>>> calcXsections(int R, TString pdfName);

//Read 2D histogram to the vector, index is rapidity, theory is given by fnlo
//Histograms have no error
vector<TH1D*> readHisto(fastNLOAlphas &fnlo)
{
    //fnlo.SetScaleFactorsMuRMuF(1.0, 1.0);
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
    //cout << "First part done" << endl;

    vector<TH1D*> hists;
    for(auto obj : bins) {
        double etaDn = obj.first;
        vector<double> binning = obj.second;

        TH1D * h = new TH1D(rn(), Form("%g", etaDn), binning.size()-1, binning.data());

        //cout << etaDn << endl;

        for(int i = 0; i < xSec.at(etaDn).size(); ++i) {
            h->SetBinContent(i+1, xSec.at(etaDn)[i]);
            h->SetBinError(i+1, 0);
        }

        hists.push_back(h);
    }

    /*
    cout << "Second part done Start" << endl;
    for(int i = 0; i < hists.size(); ++i)
        hists[i]->Print();
    cout << "Second part done End " << hists.size() <<  endl;
    */

    return hists;
}

//Get vector of histograms which includes scale unc
//(cnt, scaleUp, scaleDn)
vector<vector<TH1D*>> getScaleuncHistos(fastNLOAlphas &fnlo)
{
    fnlo.SetLHAPDFMember(0);

    vector<vector<double>> scales = {
                    {1.0, 2.0}, { 2.0, 2.0},
        {0.5, 1.0}, {1.0, 1.0}, { 2.0, 1.0},
        {0.5, 0.5}, {1.0, 0.5},
    };

    vector<vector<TH1D*>> histos;
    for(auto  s : scales) {
        fnlo.SetScaleFactorsMuRMuF(s[0], s[1]);
        auto hh = readHisto(fnlo);
        histos.push_back(readHisto(fnlo));
    }


    vector<TH1D*> hCnt(histos[0].size());
    vector<TH1D*> hUp(histos[0].size());
    vector<TH1D*> hDn(histos[0].size());

    for(int y = 0; y < histos[0].size(); ++y) { //loop over y-bins

        hCnt[y] = (TH1D*) histos[0][y]->Clone(rn());
        hUp[y]  = (TH1D*) histos[0][y]->Clone(rn());
        hDn[y]  = (TH1D*) histos[0][y]->Clone(rn());

        for(int i = 1; i <= histos[0][y]->GetNbinsX(); ++i) { //loop over pt-bins
            double cnt = histos[0][y]->GetBinContent(i);
            double up = 0, dn = 0;
            for(int s = 1; s < scales.size(); ++s) { //loop over scales
                double err  = histos[s][y]->GetBinContent(i) - cnt;
                up = max(up, err);
                dn = max(dn,-err);
            }
            hCnt[y]->SetBinContent(i, cnt);
            hUp[y]->SetBinContent(i, cnt + up);
            hDn[y]->SetBinContent(i, cnt - dn);
        }
    }

    hCnt[0]->Print();
    hUp[0]->Print();
    hDn[0]->Print();

    return {hCnt, hUp, hDn};
}

//Get histogram including up and dn pdf variation 
vector<vector<TH1D*>> getPDFuncHistos(fastNLOAlphas &fnlo)
{
    fnlo.SetLHAPDFMember(0);
    fnlo.SetScaleFactorsMuRMuF(1, 1);

    int nPDFs = fnlo.GetNPDFMembers();
    TString pdfName = fnlo.GetLHAPDFFilename();
    cout << pdfName << endl;
    double Fact = pdfName.Contains("CT14") ? 1.645 : 1;

    vector<vector<TH1D*>> histos;
    for(int i = 0; i < nPDFs; ++i) {
        fnlo.SetLHAPDFMember(i);
        histos.push_back(readHisto(fnlo));
    }

    vector<TH1D*> hCnt(histos[0].size());
    vector<TH1D*> hUp(histos[0].size());
    vector<TH1D*> hDn(histos[0].size());

    for(int y = 0; y < histos[0].size(); ++y) { //loop over y-bins

        hCnt[y] = (TH1D*) histos[0][y]->Clone(rn());
        hUp[y] = (TH1D*) histos[0][y]->Clone(rn());
        hDn[y] = (TH1D*) histos[0][y]->Clone(rn());

        for(int i = 1; i <= histos[0][y]->GetNbinsX(); ++i) { //loop over pt-bins
            double cnt = histos[0][y]->GetBinContent(i);
            double up = 0, dn = 0;
            for(int s = 1; s < nPDFs; ++s) { //loop over scales
                double err  = (histos[s][y]->GetBinContent(i) - cnt)/Fact;
                up = hypot(up, max(0.0, err));
                dn = hypot(dn, max(0.0,-err));
            }
            hCnt[y]->SetBinContent(i, cnt);
            hUp[y]->SetBinContent(i, cnt + up);
            hDn[y]->SetBinContent(i, cnt - dn);
        }
    }
    return {hCnt, hUp, hDn};
}

////rebin differential histo according to template
//TH1D *rebin(TH1D *h, TH1D *hTemp) //second is template
//{
//    TH1D *hNew = (TH1D *) hTemp->Clone(rn());
//    hNew->Reset();
//
//    for(int i = 1; i <= h->GetNbinsX(); ++i) {
//        double xCnt = h->GetBinCenter(i);
//
//        double iBin = hNew->FindBin(xCnt);
//
//        if(iBin == 0 || iBin == hNew->GetNbinsX()+1)
//            continue;
//
//        double v = h->GetBinContent(i);
//        double w = h->GetBinWidth(i);
//
//        double orgVal = hNew->GetBinContent(iBin);
//        hNew->SetBinContent(iBin, orgVal + v*w);
//    }
//
//    hNew->Scale(1, "width");
//    return hNew;
//}

//void printHisto(TH1D *h)
//{
//    for(int i = 1; i < h->GetNbinsX(); ++i) {
//        cout << i << " "<<  h->GetBinLowEdge(i) <<" "<<  h->GetXaxis()->GetBinUpEdge(i) <<" "<< h->GetBinContent(i) <<  endl;
//    }
//}

void SaveHistos(vector<vector<TH1D*>> hist,  TString tag)
{
    vector<TString> sysTag = {"Cnt", "Up", "Dn"};
    for(int s = 0; s < hist.size(); ++s) {
        for(int y = 0; y < hist[0].size(); ++y) {
            TString n = tag +"_"+ sysTag[s] +"_"+ Form("y%d", y);
            hist[s][y]->SetName(n);
            hist[s][y]->Write(n);
        }
    }
}

void setNLO(fastNLOAlphas &fnlo)
{
    fnlo.SetContributionON(fastNLO::kFixedOrder,0,true);
    fnlo.SetContributionON(fastNLO::kFixedOrder,1,true);
    fnlo.SetUnits(fastNLO::kPublicationUnits);
}

void calcXsectionsAll(int R, vector<TString> pdfNames)
{
    using namespace std;
    using namespace say;		// namespace for 'speaker.h'-verbosity levels
    using namespace fastNLO;	// namespace for fastNLO constants

    TString fastName;
    if      (R == 4) fastName = "theorFiles/InclusiveNJets_fnl5362h_v23_fix.tab"; // TODO
    else if (R == 7) fastName = "theorFiles/InclusiveNJets_fnl5332h_v23_fix.tab"; // TODO
    else exit(0);

    vector<vector<vector<TH1D*>>> histsPDF, histsScl;
    for(auto pdfName: pdfNames) {
        fastNLOAlphas fnlo(fastName.Data(), pdfName.Data(), 0);
        setNLO(fnlo);

        histsPDF.push_back(getPDFuncHistos(fnlo));
        histsScl.push_back(getScaleuncHistos(fnlo));
    }

    TFile *fOut = new TFile(Form("theorFiles/cmsJetsNLO_AK%d.root", R), "RECREATE");
    for(int i = 0; i < histsPDF.size(); ++i) {
        SaveHistos(histsPDF[i], "hist"+  pdfNames[i] + "_PDF");
        SaveHistos(histsScl[i], "hist"+  pdfNames[i] + "_Scl");
        //SaveHistos(histsAs[i],  "hist"+  pdfNames[i] + "_As");
        cout << "Saving " << pdfNames[i] << endl;
    }
    fOut->Write();
    fOut->Close();
    cout << "End " << endl;
}


//__________________________________________________________________________________________________________________________________

int main(int argc, char** argv){
	

    // namespaces
    using namespace std;
    using namespace say;		// namespace for 'speaker.h'-verbosity levels
    using namespace fastNLO;	// namespace for fastNLO constants

	SetGlobalVerbosity(ERROR);

    calcXsectionsAll(4, {"CT14nlo"/*,  "HERAPDF20_NLO_EIG",  "NNPDF31_nlo_as_0118_hessian",  "ABMP16_5_nlo",  "MMHT2014nlo68cl",
                         "CT14nnlo", "HERAPDF20_NNLO_EIG", "NNPDF31_nnlo_as_0118_hessian", "ABMP16_5_nnlo", "MMHT2014nnlo68cl"*/});

    //vector<TH1D*> readHisto(fastNLOAlphas &fnlo, TString tabName, TString pdfName)

    //say::SetGlobalVerbosity(say::DEBUG);

    //fastNLOAlphas fnloOld("theorFiles/fastnlo-cms-incjets-arxiv-1605.04436-xsec001.tab", "CT14nlo", 0);
    //fastNLOAlphas fnloNew("theorFiles/InclusiveNJets_fnl5362h_v23_fix.tab", "CT14nlo", 0);
    //fnloOld.SetContributionON(fastNLO::kFixedOrder,0,true);
    //fnloOld.SetContributionON(fastNLO::kFixedOrder,1,true);
    //fnloOld.SetUnits(fastNLO::kPublicationUnits);

    //fnloNew.SetContributionON(fastNLO::kFixedOrder,0,true);
    //fnloNew.SetContributionON(fastNLO::kFixedOrder,1,true);
    //fnloNew.SetUnits(fastNLO::kPublicationUnits);
  
    return EXIT_SUCCESS;
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
