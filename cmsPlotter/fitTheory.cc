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
#include <functional>
//#include "fastnlotk/fastNLODiffReader.h"
//#include "fastNLODiffAlphas.h"
#include "fastnlotk/fastNLOAlphas.h"

#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"


#include "plottingHelper.h"
using namespace PlottingHelper;

#include "tools.h"

// Function prototype for flexible-scale function 
double Function_Mu(double s1, double s2 );

double Function_q(double q, double pt );
double Function_pt(double q, double pt );


TString rn() {return Form("%d",rand());}

using namespace std;

const TString year = "2016";

TFile *fTh = nullptr;

vector<TH1D*> readHistos(TString pdfName, double as)
{

    vector<TH1D*> vTh;
    int idCnt = 0;
    for(int y = 0; y < 5; ++y) {
        int asI = round(as*1000);
        vTh.push_back( (TH1D*) fTh->Get(pdfName+ Form("_y%d_as0%d_scale0",y,asI) ));
    }

    for(int y = 0; y < 5; ++y) {
        applyNPEW(vTh[y],  y, year);
        applyKfactor(vTh[y], y, "kFactorNLL");
        //applyKfactor(vTh[y], y, "kFactorNNLO");
    }
    return vTh;
}

/*
vector<TH1D*> readHistos(fastNLOAlphas &fnlo)
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

    return hists;
}
*/







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


//__________________________________________________________________________________________________________________________________


struct point {
    double ptMin, ptMax;
    double yMin, yMax;
    double sigma;
    double th;
    double errStat, errSys, errTot, errUnc;
    std::vector<double> errs;//10 items
};



struct asFitter {
    vector<point> data;
    //vector<double> th;
    function<bool(point)> Cut;

    void readData(TString fName)
    {

        ifstream infile(fName);
        string line;
        bool isIn = false;
        while (getline(infile, line))
        {
            if(line.size() < 5 && line[0] == '*') {
                isIn = true;
                continue;
            }
            if(!isIn) continue;

            istringstream iss(line);
            //cout << "Line size " << line.size() << endl;
            point p;
            double flag, nPerL, nPerH, lumi;
            iss >> flag >> p.yMin >> p.yMax >> p.ptMin >> p.ptMax >> p.sigma >> p.errStat >> p.errUnc;
            p.errStat /= 100;
            p.errUnc /= 100;

            vector<double> unc;
            iss >> nPerL >> nPerH >> lumi;
            
            p.errs.push_back((nPerL-nPerH)/2);
            p.errs.push_back(lumi);
            p.th = 0;

            double a, b;
            while ((iss >> a >> b)) {
                p.errs.push_back( (a - b)/2 );
                //cout << a << " " << b << " ";
            } // error
            for(auto & e : p.errs) {
                e /= 100;
            }

            if(p.sigma > 0)
                data.push_back(p);
            //cout << endl << endl;;
            // process pair (a,b)
        }

    }


    void fillTheory(TString pdfName, double as)
    {
        vector<TH1D*> thHist = readHistos(pdfName, as);
        thHist.resize(4);

        for(auto &p : data) {
            int y = round(p.yMin * 2);
            double ptCnt = (p.ptMin + p.ptMax) / 2.;
            int binId = thHist[y]->FindBin(ptCnt);
            p.th = thHist[y]->GetBinContent(binId);
        }
        for(auto & h : thHist)
            delete h;
    }


    int getNpoints()
    {
        int s = 0;
        for(const auto &p : data)
            s += Cut(p);
        return s;
    }





    //Calculated according to https://arxiv.org/pdf/hep-ex/0012053.pdf
    //Formula (33), page 29
    double getChi2()
    {
        TVectorD s = getShifts();
        double chi2 = getChi2(s);
        
        //print shifts
        //for(int i = 0; i < data[0].errs.size(); ++i)
            //cout <<"shift " <<  i <<" "<<  s(i) << endl;
        return chi2;
    }

    //bool Cut(const point &p) { if(p.sigma == 0) return false;  return true;}


    TVectorD getShifts()
    {
        int nErr = data[0].errs.size();
        TMatrixD mat(nErr, nErr);
        TVectorD yVec(nErr);
        //Calculate the optimal shifts
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double th = p.th;
            double C = pow(p.sigma*p.errStat,2) + pow(p.sigma*p.errUnc,2);
            for(int j = 0; j < p.errs.size(); ++j)
                for(int k = 0; k < p.errs.size(); ++k)
                    mat(j,k) += 1./C * th*th * p.errs[j]*p.errs[k];

            for(int j = 0; j < p.errs.size(); ++j)
                yVec(j) += - 1./C * (p.sigma - th) * th * p.errs[j];

        }

        for(int j = 0; j < data[0].errs.size(); ++j)
            mat(j,j) += 1;

        //Solve 
        TDecompSVD svd(mat);
        Bool_t ok;
        const TVectorD sh = svd.Solve(yVec, ok);

        return sh;
    }

    double getChi2(const TVectorD &s)
    {
        //Evaluate the chi2
        double chi2 = 0;
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double th = p.th;
            double corErr = 0;
            for(int i = 0; i < p.errs.size(); ++i)
                corErr += s(i) * p.errs[i];

            double C = pow(p.sigma*p.errStat,2) + pow(p.sigma*p.errUnc,2);
            chi2 += pow(p.sigma - th * (1 - corErr), 2) / C;
        }

        for(int j = 0; j < data[0].errs.size(); ++j)
            chi2 += pow(s(j),2);

        return chi2;
    }

    /*
    //Fill theory
    vector<TH1D*> readHistos(fastNLOAlphas &fnlo)
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


        //cout << "Second part done Start" << endl;
        //for(int i = 0; i < hists.size(); ++i)
        //    hists[i]->Print();
        //cout << "Second part done End " << hists.size() <<  endl;


        return hists;
    }
    */

    double calcChi2(TString pdfName, double as) {

        fillTheory(pdfName, as);

        return getChi2();
    }

    void printChi2Table()
    {
        vector<TString> pdfSets = {"CT14nnlo", "HERAPDF20_NNLO", "NNPDF31_nnlo"};

        for(auto p : pdfSets)
            cout << p << " ";
        cout << endl;

        for(int y = 0; y < 4; ++y) { //over rapidities
            cout << y*0.5 <<" " << (y+1)*0.5 << " ";
            Cut = [y](point p) { return ( abs(y*0.5-p.yMin) < 0.1 &&  p.sigma != 0);};
            int ndf = getNpoints();

            for(auto pdfSet : pdfSets) {
                double chi2now = calcChi2(pdfSet, 0.118);
                cout  <<chi2now << " / " << ndf << " " ;
            }
            cout << endl;

        }
        //double chi2now = asfit.calcChi2("CT14nnlo", 0.118);
        //cout << as <<" : "<<chi2now << " / " << ndf << endl;
    }

    void printAs(int y)
    {
        Cut = [y](point p) { return ( abs(y*0.5-p.yMin) < 0.1 &&  p.sigma != 0);};
        int ndf = getNpoints();
        for(double as = 0.113; as <= 0.122; as +=0.001) {
            double chi2now = calcChi2("HERAPDF20_NNLO", as);
            cout << y <<" "<< as <<" : "<<chi2now << " / " << ndf << endl;
        }
    }




};





int main(int argc, char** argv)
{
    fTh  = TFile::Open("cmsJetsAsScan.root");  //NLO predictions

	asFitter asfit;
    asfit.readData("xFitterTables/data16.txt");


    asfit.printAs(0);
    asfit.printAs(1);
    asfit.printAs(2);

    //asfit.printChi2Table();
    return 0;


    for(double as = 0.113; as <= 0.122; as +=0.001) {
        asfit.Cut = [](point p) { return (p.sigma != 0);};
        int ndf = asfit.getNpoints();
        double chi2now = asfit.calcChi2("CT14nnlo", 0.118);
        cout << as <<" : "<<chi2now << " / " << ndf << endl;
    }

    return 0;

  
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
