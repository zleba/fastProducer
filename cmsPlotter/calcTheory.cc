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

//Read 2D histogram to the vector, index is rapidity, theory is given by fnlo
//Histograms have no error
TH1D* readHisto(fastNLOAlphas &fnlo)
{
    //fnlo.SetScaleFactorsMuRMuF(1.0, 1.0); // this is done before calling `readHisto`
    fnlo.CalcCrossSection();
    vector<double> xs = fnlo.GetCrossSection();
    
    vector<double> bins;
    for(int k = 0; k < xs.size(); ++k) {
        //cout << "Helenka " << k <<" "<< fnlodiff.GetObsBinLoBound(k,0) << endl;
        double ptDn  = fnlo.GetObsBinLoBound(k,0);
        double ptUp  = fnlo.GetObsBinUpBound(k,0);

        //cout << setw(10) << ptDn << setw(10) << ptUp << setw(15) << xs.at(k) << '\n';

        // bin edges
        bins.push_back(ptDn);
        if (k == xs.size() - 1 || fnlo.GetObsBinLoBound(k,0) != fnlo.GetObsBinLoBound(k+1,0))
            bins.push_back(ptUp);
    }
    //cout << flush;

    TH1D * h = new TH1D(rn(), "", bins.size()-1, bins.data());

    for (int i = 0; i < xs.size(); ++i) {
        h->SetBinContent(i+1, xs[i]);
        h->SetBinError(i+1, 0);
    }

    return h;
}

//Get vector of histograms which includes scale unc
//(cnt, scaleUp, scaleDn)
vector<TH1D*> getScaleuncHistos(fastNLOAlphas &fnlo)
{
    fnlo.SetLHAPDFMember(0);

    vector<vector<double>> scales = {
                    {1.0, 2.0}, { 2.0, 2.0},
        {0.5, 1.0}, {1.0, 1.0}, { 2.0, 1.0},
        {0.5, 0.5}, {1.0, 0.5},
    };

    vector<TH1D*> histos;
    for(auto  s : scales) {
        fnlo.SetScaleFactorsMuRMuF(s[0], s[1]);
        auto hh = readHisto(fnlo);
        histos.push_back(readHisto(fnlo));
    }


    auto hCnt = (TH1D*) histos[0]->Clone(rn());
    auto hUp  = (TH1D*) histos[0]->Clone(rn());
    auto hDn  = (TH1D*) histos[0]->Clone(rn());

    for(int i = 1; i <= histos[0]->GetNbinsX(); ++i) { //loop over pt-bins
        double cnt = histos[0]->GetBinContent(i);
        double up = 0, dn = 0;
        for(int s = 1; s < scales.size(); ++s) { //loop over scales
            double err  = histos[s]->GetBinContent(i) - cnt;
            up = max(up, err);
            dn = max(dn,-err);
        }
        hCnt->SetBinContent(i, cnt);
        hUp->SetBinContent(i, cnt + up);
        hDn->SetBinContent(i, cnt - dn);
    }

    hCnt->Print();
    hUp->Print();
    hDn->Print();

    return {hCnt, hUp, hDn};
}

//Get histogram including up and dn pdf variation 
vector<TH1D*> getPDFuncHistos(fastNLOAlphas &fnlo)
{
    fnlo.SetLHAPDFMember(0);
    fnlo.SetScaleFactorsMuRMuF(1, 1);

    int nPDFs = fnlo.GetNPDFMembers();
    TString pdfName = fnlo.GetLHAPDFFilename();
    double Fact = pdfName.Contains("CT14") ? 1.645 : 1; // TODO?

    vector<TH1D*> histos;
    for(int i = 0; i < nPDFs; ++i) {
        fnlo.SetLHAPDFMember(i);
        histos.push_back(readHisto(fnlo));
    }

    auto hCnt = (TH1D*) histos[0]->Clone(rn());
    auto hUp = (TH1D*) histos[0]->Clone(rn());
    auto hDn = (TH1D*) histos[0]->Clone(rn());

    for(int i = 1; i <= histos[0]->GetNbinsX(); ++i) { //loop over pt-bins
        double cnt = histos[0]->GetBinContent(i);
        double up = 0, dn = 0;
        for(int s = 1; s < nPDFs; ++s) { //loop over scales
            double err  = (histos[s]->GetBinContent(i) - cnt)/Fact;
            up = hypot(up, max(0.0, err));
            dn = hypot(dn, max(0.0,-err));
        }
        hCnt->SetBinContent(i, cnt);
        hUp->SetBinContent(i, cnt + up);
        hDn->SetBinContent(i, cnt - dn);
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

    TFile *fOut = new TFile(Form("theorFiles/cmsJetsNLO_AK%d.root", R), "RECREATE");

    for (int y = 0; y < 5; ++y) { // TODO: 5
        fOut->cd();
        auto yDir = fOut->mkdir(Form("ybin%d", y+1));
        yDir->cd();

        cout << "y = " << y << endl;

        TString fastName = Form("data/ak%d/NNLO/1jet.NNLO.fnl5362h_y%d_ptjet.tab.gz", R, y);
        cout << fastName << endl;
        //if      (R == 4) fastName = "theorFiles/InclusiveNJets_fnl5362h_v23_fix.tab";
        //else if (R == 7) fastName = "theorFiles/InclusiveNJets_fnl5332h_v23_fix.tab";
        //else exit(0);

        vector<vector<TH1D*>> histsPDF, histsScl;
        for(auto pdfName: pdfNames) {

            yDir->cd();
            auto pdfDir = yDir->mkdir(pdfName);
            pdfDir->cd();

            cout << pdfName << endl;
            fastNLOAlphas fnlo(fastName.Data(), pdfName.Data(), 0);
            fnlo.SetUnits(fastNLO::kPublicationUnits);
            fnlo.SetContributionON(fastNLO::kFixedOrder, 0, true); // LO
            fnlo.SetContributionON(fastNLO::kFixedOrder, 1, true); // NLO
            bool doNNLO = pdfName.Contains("NNLO") || pdfName.Contains("nnlo");
            fnlo.SetContributionON(fastNLO::kFixedOrder, 2, doNNLO); // NNLO

            auto histPDF = getPDFuncHistos  (fnlo);
            auto histScl = getScaleuncHistos(fnlo);

            histPDF.at(0)->Write("nominal");
            histPDF.at(1)->Write("PDFup");
            histPDF.at(2)->Write("PDFdn");
            histScl.at(1)->Write("Sclup");
            histScl.at(2)->Write("Scldn");
        }
    }
    fOut->Close();
    cout << "End " << endl;
}


//__________________________________________________________________________________________________________________________________

int main(int argc, char** argv)
{

    // namespaces
    using namespace std;
    using namespace say;		// namespace for 'speaker.h'-verbosity levels
    using namespace fastNLO;	// namespace for fastNLO constants

	SetGlobalVerbosity(ERROR);

    calcXsectionsAll(4, {"CT14nlo",  "HERAPDF20_NLO_EIG",  "NNPDF31_nlo_as_0118_hessian",  "ABMP16_5_nlo",  "MMHT2014nlo68cl",
                         "CT14nnlo", "HERAPDF20_NNLO_EIG", "NNPDF31_nnlo_as_0118_hessian", "ABMP16_5_nnlo", "MMHT2014nnlo68cl"});

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
