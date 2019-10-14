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
#include "Math/Functions.h"
#include "TF1.h"


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
    std::vector<double> errs, thErrs;//10 items
};



struct asFitter {
    vector<point> data; //allDataPoints + potential theory predictions
    //vector<double> th;
    function<bool(point)> Cut; //selection function

    //Map with theorXsections [pdfName][alphaS*1000] [scaleVar][iPdf][rap]
    map<TString, map<int, vector< vector<vector<TH1D*>> >>>  thHists; 

    //Read data from the text file
    static vector<point>  readData(TString fName)
    {
        vector<point> dataNow;

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
                dataNow.push_back(p);
            //cout << endl << endl;;
            // process pair (a,b)
        }
        return dataNow;

    }


    //Read theory histogram pdfName, as and scale variation s (the NP/EW corrections are applied)
    static vector<vector<TH1D*>> readHistos(TString pdfName, double as, int s)
    {

        //Size of the histogram depends on the as value and #pdf for given pdf set
        
        vector<vector<TH1D*>> vTh; //indexes - [pdfVar][y]
        vTh.resize(1);
        if(pdfName.Contains("CT14")  && abs(as - 0.118) < 1e-6)
            vTh.resize(56+1);
        else if(pdfName.Contains("HERAPDF")  && abs(as - 0.118) < 1e-6)
            vTh.resize(28+1);
        else if(pdfName.Contains("NNPDF31")  && abs(as - 0.118) < 1e-6)
            vTh.resize(100+1);

        int idCnt = 0;
        for(int ipdf = 0; ipdf < vTh.size(); ++ipdf) {
            for(int y = 0; y < 5; ++y) {
                int asI = round(as*1000);
                vTh[ipdf].push_back( (TH1D*) fTh->Get(pdfName+ Form("_y%d_as0%d_scale%d_pdf%d",y,asI, s, ipdf) )->Clone(rn())  );
            }
        }

        for(int ipdf = 0; ipdf < vTh.size(); ++ipdf) {
            for(int y = 0; y < 5; ++y) {
                applyNPEW(vTh[ipdf][y],  y, year);
                applyKfactor(vTh[ipdf][y], y, "kFactorNLL");
                //applyKfactor(vTh[y], y, "kFactorNNLO");
            }
        }

        return vTh;
    }

    //Read theory histograms for PDF pdfName (all alphaS (as) and all scale choices (s))
    void readAllTheory(TString pdfName) {
        for(auto as: pdfAsVals.at(pdfName)) {
            cout << pdfName <<" "<< as << endl;
            int asI = round(as*1000);
            thHists[pdfName][asI].resize(7);
            for(int s = 0; s < 7; ++s)
                thHists[pdfName][asI][s] = readHistos(pdfName, as, s);
        }
    }


    //Fill theory to the points in vector<points>, resutl contains also PDF unc.
    void fillTheory(TString pdfName, double as, int scale = 0)
    {
        //vector<vector<TH1D*>> thHist    = readHistos(pdfName, as);
        //vector<vector<TH1D*>> thHist118 = readHistos(pdfName, 0.118);
        
        //cout << pdfName <<" "<< as <<" : begin"<< endl;
        int asI = round(as*1000);
        vector<vector<TH1D*>> thHist    = thHists.at(pdfName).at(asI)[scale];
        vector<vector<TH1D*>> thHist118 = thHists.at(pdfName).at(118)[scale];
        //cout << pdfName <<" "<< as <<" : end"<< endl;

        //thHist.resize(4);

        for(auto &p : data) {
            int y = round(p.yMin * 2);
            double ptCnt = (p.ptMin + p.ptMax) / 2.;
            int binId = thHist[0][y]->FindBin(ptCnt);
            p.th = thHist[0][y]->GetBinContent(binId);

            p.thErrs.clear();

            for(int i = 1; i < thHist118.size(); ++i) {
                double thNom = thHist118[0][y]->GetBinContent(binId);
                double diff = (thHist118[i][y]->GetBinContent(binId) - thNom) / thNom;
                p.thErrs.push_back(diff);
            }


        }
    }


    //Get number of points fulfilling the cuts
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


    //Get the vector with the nuissence parameters (values which minimize chi2)
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

    //get the chi2 vale, the nuisence vector is as an input
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

    //Get chi2 based on covariance matrix
    double getChi2cov()
    {
        //Filter data
        vector<point> dataF;
        for(const auto &p : data) {
            if(!Cut(p)) continue;
            dataF.push_back(p);
        }

        //Fill stat cov matrix
        TMatrixD CovStat(dataF.size(), dataF.size());
        for(int i = 0; i < dataF.size(); ++i) {
            const auto &p = dataF[i];
            double C = pow(p.sigma*p.errStat,2) + pow(p.sigma*p.errUnc,2);
            CovStat(i,i) = C;
        }

        //Fill data sys matrix
        TMatrixD CovSys(dataF.size(), dataF.size());
        for(int k = 0; k < dataF[0].errs.size(); ++k) {
            for(int i = 0; i < dataF.size(); ++i) 
            for(int j = 0; j < dataF.size(); ++j) {
                //if(dataF[i].yMin != dataF[j].yMin)// && k <= 15)
                 //   continue;
                CovSys(i,j) += dataF[i].errs[k]*dataF[j].errs[k]  *  dataF[i].sigma * dataF[j].sigma;
            }
        }

        //Fill pdf sys matrix

        double Ccorr = 1;
        if(dataF[0].thErrs.size() == 57) //CT14
            Ccorr = 1*1./(2*1.64*1.64);
        else if(dataF[0].thErrs.size() == 29) //HERAPDF
            Ccorr = 1*1./(2);
        else if(dataF[0].thErrs.size() == 101) //NNPDF31
            Ccorr = 1*1./(2);

        TMatrixD CovPDF(dataF.size(), dataF.size());
        for(int k = 0; k < dataF[0].thErrs.size(); ++k) {
            for(int i = 0; i < dataF.size(); ++i) 
            for(int j = 0; j < dataF.size(); ++j) {
                CovPDF(i,j) += dataF[i].thErrs[k]*dataF[j].thErrs[k]  *  dataF[i].sigma * dataF[j].sigma * Ccorr;
            }
        }

        TMatrixD Cov = CovStat + CovSys +  CovPDF;// + CovSys + CovPDF;

        double chi2 = 0;
        
        //Fill data - th
        TVectorD diff(dataF.size());
        for(int i = 0; i < dataF.size(); ++i) {
            const auto &p = dataF[i];
            diff(i) = p.sigma - p.th;
        }


        //Evaluate chi2
        TMatrixD CovInv = Cov;
        CovInv.Invert();
        TVectorD diffM = CovInv*diff;

        for(int i = 0; i < dataF.size(); ++i)
            chi2 += diff(i) * diffM(i);


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

    //fill the theory from file to histos and get chi2 wrt data
    double calcChi2(TString pdfName, double as, int scale = 0) {

        fillTheory(pdfName, as, scale);

        return getChi2cov();
    }


    //print the chi2 table
    void printChi2Table()
    {
        vector<TString> pdfSets = {"CT14nlo", "HERAPDF20_NLO", "NNPDF31_nnlo"};

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

        
        { //Total chi2 
            cout << "Total ";
            Cut = [](point p) { return ( abs(p.yMin) < 1.6 &&  p.sigma != 0);};
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

    void printAsY(int y)
    {
        if(y < 0) Cut = [y](point p) { return ( abs(p.yMin) < 1.6 &&  p.sigma != 0);};
        else      Cut = [y](point p) { return ( abs(y*0.5-p.yMin) < 0.1 &&  p.sigma != 0);};

        int ndf = getNpoints();
        for(double as = 0.113; as <= 0.122; as +=0.001) {
            //double chi2now = calcChi2("HERAPDF20_NLO", as);
            double chi2now = calcChi2("CT14nlo", as);
            cout << y <<" "<< as <<" : "<<chi2now << " / " << ndf << endl;
        }
    }

    void printAsPt(int pt)
    {
        const vector<double> ptBinsAs = {97, 174, 272, 395, 548, 737, 967, 1248, 1588, 2000, 2500, 3103};
        double ptMin = ptBinsAs[pt];
        double ptMax = ptBinsAs[pt+1];
        Cut = [ptMin,ptMax](point p) { return (  abs(p.yMin) < 1.6 &&  p.sigma != 0   && p.ptMin >= ptMin -1 &&  p.ptMax <= ptMax +1     );};


        int ndf = getNpoints();
        for(double as = 0.113; as <= 0.122; as +=0.001) {
            //double chi2now = calcChi2("HERAPDF20_NLO", as);
            double chi2now = calcChi2("CT14nlo", as);
            cout << ptBinsAs[pt] <<" "<< as <<" : "<<chi2now << " / " << ndf << endl;
        }
    }

    TGraph *getFitGraphPt(TString pdfName, int pt, int scale)
    {
        double ptMin = ptBinsAs[pt];
        double ptMax = ptBinsAs[pt+1];
        Cut = [ptMin,ptMax](point p) { return (  abs(p.yMin) < 0.3 &&  p.sigma != 0   && p.ptMin >= ptMin -1 &&  p.ptMax <= ptMax +1     );};

        int ndf = getNpoints();

        TGraph *gr = new TGraph();

        int i = 0;
        for(double as  : pdfAsVals.at(pdfName) ) {
            double chi2now = calcChi2(pdfName, as, scale);
            //cout << y <<" "<< as <<", scale="<<scale <<" : "<<chi2now << " / " << ndf << endl;
            gr->SetPoint(i, as, chi2now);
            ++i;
        }
        gr->Fit("pol4");
        return gr;
    }




    TGraph *getFitGraph(TString pdfName, int y, int scale)
    {
        if(y < 0) Cut = [y](point p) { return ( abs(p.yMin) < 1.6 &&  p.sigma != 0);};
        else      Cut = [y](point p) { return ( abs(y*0.5-p.yMin) < 0.1 &&  p.sigma != 0);};
        //Cut = [y](point p) { return ( abs(p.yMin) < 1.6 &&  p.sigma != 0);};
        int ndf = getNpoints();

        TGraph *gr = new TGraph();

        int i = 0;
        for(double as  : pdfAsVals.at(pdfName) ) {
            double chi2now = calcChi2(pdfName, as, scale);
            cout << y <<" "<< as <<", scale="<<scale <<" : "<<chi2now << " / " << ndf << endl;
            gr->SetPoint(i, as, chi2now);
            ++i;
        }
        gr->Fit("pol4");
        return gr;
    }


    void getAllChi2s()
    {

        vector<TString> pdfNames;
        for(auto el :  thHists)
            pdfNames.push_back(el.first);

        TFile *fOut = TFile::Open("chi2.root", "RECREATE");
        for(auto pdfName : pdfNames) {
            for(int y = -1; y < 4; ++y) {
                for(int s = 0; s < 7; ++s) {
                    TGraph *gr = getFitGraph(pdfName, y, s);
                    if(y == -1) gr->Write(pdfName + Form("_scale%d", s));
                    else gr->Write(pdfName + Form("_Y%d_scale%d", y, s));
                }
            }
        }

        //Fill pT dep
        for(auto pdfName : pdfNames) {
            for(int pt = 0; pt < ptBinsAs.size()-1; ++pt) {
                for(int s = 0; s < 7; ++s) {
                    TGraph *gr = getFitGraphPt(pdfName, pt, s);
                    gr->Write(pdfName + Form("_Pt%d_scale%d", pt, s));
                }
            }
        }


        fOut->Write();
        fOut->Close();
    }



    void fitAs(TString pdfName, int y)
    {
        for(int s = 0; s < 7; ++s) {
            TCanvas *c = new TCanvas(rn(), "can", 600, 600);
            TGraph *gr = getFitGraph(pdfName, y, s);

            TF1 *fit = gr->GetFunction("pol4");
            double asMin = fit->GetMinimumX();
            double chi2Min = fit->Eval(asMin);
            double shLow  = fit->GetX(chi2Min + 1, 0, asMin);
            double shHigh = fit->GetX(chi2Min + 1, asMin, 1);

            double errL = asMin - shLow;
            double errH = shHigh - asMin;
            cout << "Helenka min " << asMin << " "<< errL <<" "<< errH << endl;
            gr->Draw("a*");
            c->SaveAs(Form("asFit%d.pdf", s));
        }

    }


};





int main(int argc, char** argv)
{
    fTh  = TFile::Open("cmsJetsAsScan.root");  //NLO predictions

	asFitter asfit;
    asfit.data = asfit.readData("xFitterTables/data16.txt");
    asfit.readAllTheory("CT14nlo");
    asfit.readAllTheory("HERAPDF20_NLO");
    asfit.readAllTheory("NNPDF31_nnlo");


    asfit.printChi2Table();

    return 0;
    asfit.getAllChi2s();
    return 0;

    asfit.Cut = [](point p) { return ( abs(p.yMin) < 1.7 &&  p.sigma != 0);};
    for(auto as: pdfAsVals.at("CT14nlo")) {
        cout << "Helenka " << as <<" "<<  asfit.calcChi2("CT14nlo", as, 0) <<" : "<< asfit.getNpoints() << endl;
    }

    return 0;

    for(int y = 0; y < 1; ++y) {
        asfit.Cut = [y](point p) { return ( abs(y*0.5-p.yMin) < 0.1 &&  p.sigma != 0);};
        for(auto as: pdfAsVals.at("CT14nlo")) {
            cout << "Helenka " << y <<" "<<as <<" "<<  asfit.calcChi2("CT14nlo", as, 0) <<" : "<< asfit.getNpoints() << endl;
        }
    }




    return 0;


    asfit.fitAs("CT14nlo", -1);
    return 0;

    for(int i = 0; i < 11; ++i)
        asfit.printAsPt(i);

    return 0;


    asfit.printAsY(-1);
    asfit.printAsY(0);
    asfit.printAsY(1);
    asfit.printAsY(2);
    asfit.printAsY(3);

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
