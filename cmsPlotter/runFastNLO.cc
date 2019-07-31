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

        for(int i = 0; i < xSec.at(etaDn).size(); ++i)
            h->SetBinContent(i+1, xSec.at(etaDn)[i]);

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

int main(int argc, char** argv){
	

  // namespaces
  using namespace std;
  using namespace say;		// namespace for 'speaker.h'-verbosity levels
  using namespace fastNLO;	// namespace for fastNLO constants

	SetGlobalVerbosity(DEBUG);


    //vector<TH1D*> readHisto(fastNLOLHAPDF &fnlo, TString tabName, TString pdfName)

    //say::SetGlobalVerbosity(say::DEBUG);
    fastNLOLHAPDF fnloOld("fastTables/fastnlo-cms-incjets-arxiv-1605.04436-xsec001.tab", "CT10", 0);
    fastNLOLHAPDF fnloNew("fastTables/InclusiveNJets_fnl5362h_v23_fix.tab", "CT10", 0);
    fnloOld.SetContributionON(fastNLO::kFixedOrder,0,true);
    fnloOld.SetContributionON(fastNLO::kFixedOrder,1,true);
    fnloOld.SetUnits(fastNLO::kPublicationUnits);

    fnloNew.SetContributionON(fastNLO::kFixedOrder,0,true);
    fnloNew.SetContributionON(fastNLO::kFixedOrder,1,true);
    fnloNew.SetUnits(fastNLO::kPublicationUnits);


    int nMem = fnloOld.GetNPDFMembers();


    vector<TH1D*> oldH = readHisto(fnloOld);
    vector<TH1D*> newH = readHisto(fnloNew);



    int yB = 2;
    TH1D *oldHH = rebin(oldH[yB], newH[yB]);

    for(int i = 1; i < oldHH->GetNbinsX(); ++i) {
        cout << i << " "<<  oldHH->GetBinLowEdge(i) <<" "<<  oldHH->GetBinContent(i) <<" "<< newH[yB]->GetBinContent(i) <<  endl;
    }

    return 0;

    TCanvas *can = new TCanvas("can", "");
    can->SetLogx();
    can->SetLogy();
    
    oldH[0]->Draw();
    newH[0]->Draw("same");

    can->SaveAs("plot.pdf");



    /*

    vector<vector<double>> xs(9 + nMem); //scale variations + PDF variations

    double vars[] = {1, 2., 0.5};

    //Scale variations 
    for(int iR = 0; iR <= 2; ++iR)
    for(int iF = 0; iF <= 2; ++iF) {
        int iG = 3*iR + iF;
        fnlodiff.SetScaleFactorsMuRMuF(vars[iR], vars[iF]);
        fnlodiff.CalcCrossSection();
        xs[iG] = fnlodiff.GetCrossSection();
    }

    //Scale variations 



    for(int k = 0; k < xs[0].size(); ++k) {
        //cout << "Helenka " << k <<" "<< fnlodiff.GetObsBinLoBound(k,0) << endl;
        double etaAvg = (fnlodiff.GetObsBinLoBound(k,0) + fnlodiff.GetObsBinUpBound(k,0))/2;
        double ptAvg = (fnlodiff.GetObsBinLoBound(k,1) + fnlodiff.GetObsBinUpBound(k,1))/2;

        double w = fnlodiff.GetObsBinUpBound(k,0) - fnlodiff.GetObsBinLoBound(k,0);
        cout << etaAvg <<" "<<ptAvg << " "<< xs[0][k] <<" "<< xs[1][k] <<" "<< xs[2][k] <<" "<< xs[2][k]<<  endl;
    }

    */


  /*

  //  If you want to receive your cross section in
  //   pb/GeV or in pb. Here we choose pb/GeV
  fnlodiff.SetUnits(fastNLO::kPublicationUnits);
   fnlodiff.SetFit(fastNLODiffAlphas::ABCDE);
    fnlodiff.SettIntegratedRange(-1.);
    fnlodiff.SetProtonE(920.);
  
	fnlodiff.SetContributionON(fastNLO::kFixedOrder,0,true);
	fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,true);
	//fnlodiff.SetContributionON(fastNLO::kFixedOrder,2,false);
	fnlodiff.SetContributionON(fastNLO::kFixedOrder,2,true);

    
	cout << "The file name is " << argv[1] << endl;


  // Set the xpom integration interval and method
  // -------- Boris LRG dijets
  //fnlodiff.SetXPomLinSlicing( 5, xpomBins[i],  xpomBins[i+1]); // e.g. Boris LRG dijets
  fnlodiff.SetXPomLinSlicing( 30, 0.000 ,  0.030 ); // e.g. Radek VFPS dijets

  fnlodiff.SetExternalFuncForMuR (&Function_Mu);
  fnlodiff.SetExternalFuncForMuF (&Function_Mu);
  fnlodiff.SettIntegratedRange(-1.);


  vector<double> bins = { 
      0.0032,
      0.0063,
      0.0126,
      0.0251,
      0.0501,
      0.1000,
      0.1995,
      0.3981};

  TH1D *hBeta = new TH1D("hBeta", "Beta", bins.size()-1, bins.data());

  int nBins = 120;
  double step = 0.03 / nBins;


  for(int i = 0; i < nBins; ++i) {
    double xpLow = i*step;
    double xpHi  = (i+1)*step;
    
    double xpavg = (xpLow + xpHi) / 2;
   
    fnlodiff.SetXPomLinSlicing( 1, xpLow,  xpHi); // e.g. Boris LRG dijets

    vector<double>  xs = fnlodiff.GetDiffCrossSection();
    for(int k = 0; k < xs.size(); ++k) {
        cout << "Helenka " << k <<" "<< fnlodiff.GetObsBinLoBound(k,0) << endl;
        double xBj = (fnlodiff.GetObsBinLoBound(k,0) + fnlodiff.GetObsBinUpBound(k,0))/2;

        double w = fnlodiff.GetObsBinUpBound(k,0) - fnlodiff.GetObsBinLoBound(k,0);

        double beta = xBj / xpavg;

        hBeta->Fill(beta, w*xs[k]);
    }
  }

  hBeta->Scale(1.0, "width");

  for(int i = 1; i <= hBeta->GetNbinsX(); ++i)
     cout << i <<" "<< hBeta->GetBinContent(i) << endl;
  
  return 0;



  //fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,false);
  
    //fnlodiff.SetExternalFuncForMuR (&Function_Mu);
    //fnlodiff.SetExternalFuncForMuF (&Function_Mu);


  //   fnlodiff.SetAlphasMz(0.1180);

  //  const int nStep = 1;
  

//    double xpom[nStep] = {0.0299028233407105 };
//    double dxpom[nStep] = {0.000315879371143547};

 
//  fnlodiff.SetXPomSlicing( nStep, xpom, dxpom );

  /*
  // -------- Radek VFPS dijets
  fnlodiff.SetExternalFuncForMuR (&Function_Mu);
  fnlodiff.SetExternalFuncForMuF (&Function_Mu);
  fnlodiff.SettIntegratedRange(-0.6);
  fnlodiff.SetXPomLinSlicing( 30, 0.010 ,  0.024 ); // e.g. Radek VFPS dijets
  
  // Radek's slicing
//   int nStep = 56;
//   double xpom[56] = { 0.01025, 0.01075, 0.01125, 0.01175, 0.01225, 0.01275, 0.01325, 0.01375, 0.01425, 0.01475, 0.01525, 0.01575, 0.01625, 0.01675, 0.01725, 0.01775, 0.01825, 0.01875, 0.01925, 0.01975, 0.02025, 0.02075, 0.02125, 0.02175, 0.02225, 0.02275,  0.02325, 0.02375,  };
//   double dxpom[] = {0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005 , 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005};
//   fnlodiff.SetXPomSlicing( nStep, xpom, dxpom );
  
  // switch NLO contribution off
  //fnlodiff.SetContributionON(fastNLO::kFixedOrder,1,false);
     
  // make your scale definition (see above)
  //fnlodiff.SetScaleFactorsMuRMuF(1.0,1.0);
  
  // calculate and access the cross section
  vector<double>  xs = fnlodiff.GetDiffCrossSection();
  fnlodiff.PrintCrossSections();
  cout << "Print out finished" << endl;

  return 0;
  //fnlodiff.PrintCrossSectionsWithReference();

  cout<<"double[5] = {"
      <<xs[0]<<", "
      <<xs[1]<<", "
      <<xs[2]<<", "
      <<xs[3]<<", "
      <<xs[4]<<"};"<<endl;
     

  double xSec = 0;
  for(unsigned i = 0; i < xs.size(); ++i) {
      double s = abs( fnlodiff.GetObsBinLoBound(i,0) - fnlodiff.GetObsBinUpBound(i,0) );
      xSec += s * xs[i];
  }
  //cout << "Total is " << xpomBins[i]<<" "<<  xSec / (xpomBins[i+1] - xpomBins[i]) * hadCorr[i] / 1.2  << endl;

  //}


  return 0;


  
  const int nbins = fnlodiff.GetNObsBin();
  const int npdfall = 30;
  // const int npdf = 30;
  vector<vector<double> > xsPDFerr;
  xsPDFerr.push_back(xs);
  for ( int i = 1; i<=npdfall ; i++ ) {
     fnlodiff.SetLHAPDFMember(i);
     //fnlodiff.CalcCrossSection();
     vector<double>  xsErr = fnlodiff.GetDiffCrossSection();
     xsPDFerr.push_back(xsErr);
     for ( int b = 0 ; b<nbins ; b++ ) {
	cout<<"  "<<xsErr[b]/xs[b]*100.<<endl;
     }
     cout<<endl;
  }

  
  vector<double> PDFerrUp, PDFerrDn, PDFerrMasterUp, PDFerrMasterDn;
//   cout<<" bin     \t\tDPDF error"<<endl;
//   for ( int i = 1; i<=npdf ; i++ ) {
//      cout<<"\n-------- PDFset "<<i<<" --------"<<endl;
//       for ( int b = 0 ; b<nbins ; b++ ) {
// 	 //cout<<"  "<<b<<"\t"<<xsPDFerr[i][b]<<"\t\t"<<xsPDFerr[i][b]/xsPDFerr[0][b]<<endl;
// 	 cout<<"  "<<b<<"\t\t"<<xsPDFerr[i][b]/xsPDFerr[0][b]<<endl;
//       }    
//   }  
  //   cout<<"\n===================================="<<endl<<endl;
  //   for ( int i = npdf+1; i<=npdfall ; i++ ) {
  //      cout<<"\n-------- PDFset "<<i<<" --------"<<endl;
  //       for ( int b = 0 ; b<nbins ; b++ ) {
  // 	 cout<<"  "<<b<<"\t"<<xsPDFerr[i][b]<<"\t\t"<<xsPDFerr[i][b]/xsPDFerr[0][b]<<endl;
  //       }    
  //   }  

  //zpom 0.2 0.4 646.5 646.1 646.3 643.6 648.4 643.3 648.2 663.8 628.1 641.4 622 603 659 644.5 647.5 682 612.2 613.6 676 643.9 647.4 643.6 648 583.8 667.8 653.7 637.1 650.6 640.1 610.7 681 
  cout<<endl;
  for ( unsigned int b = 0 ; b<xs.size() ; b++ ) {
      cout<<argv[1];
      //     cout<<" "<<fnlodiff.GetLoBin(b,0)<<" "<<fnlodiff.GetUpBin(b,0);
      cout<<" "<<fnlodiff.GetObsBinLoBound(b,0)<<" "<<fnlodiff.GetObsBinUpBound(b,0);
      for ( int i = 0; i<=npdfall ; i++ ) {
          cout<<" "<<xsPDFerr[i][b];
      }
      cout<<endl;
  }
  cout<<endl;
  */
  
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
