//********************************************************************
//     
//     fnlo-tk-h1diffpdf.cc
//     Program to read fastNLO v2 tables and derive
//     QCD cross sections using PDFs e.g. from LHAPDF
//     
//********************************************************************
#include <iostream>
#include <vector>
#include <set>
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

/*
const vector<TString> ErrNames = {
"nperr", "lumi", "AbsoluteStat",  "AbsoluteScale",  "AbsoluteMPFBias",  "Fragmentation",  "SinglePionECAL",  "SinglePionHCAL",  "FlavorQCD",  "TimePtEta",  "RelativeJEREC1",  "RelativeJEREC2",  "RelativeJERHF",  "RelativePtBB",  "RelativePtEC1",  "RelativePtEC2", "RelativePtHF",  "RelativeBal",  "RelativeSample",  "RelativeFSR",  "RelativeStatFSR",  "RelativeStatEC",  "RelativeStatHF",  "PileUpDataMC",  "PileUpPtRef",  "PileUpPtBB",  "PileUpPtEC1",  "PileUpPtEC2",  "PileUpPtHF",  "fake",  "miss",  "JER",  "PUprof"};
*/

const vector<TString> ErrNames = {
"NPerr", "Lumi", "AbsStat",  "AbsScale",  "AbsMPFBias",  "Frag",  "SinglePionECAL",  "SinglePionHCAL",  "FlavorQCD",  "TimePtEta",  "RelJEREC1",  "RelJEREC2",  "RelJERHF",  "RelPtBB",  "RelPtEC1",  "RelPtEC2", "RelPtHF",  "RelBal",  "RelSample",  "RelFSR",  "RelStatFSR",  "RelStatEC",  "RelStatHF",  "PUDataMC",  "PUPtRef",  "PUPtBB",  "PUPtEC1",  "PUPtEC2",  "PUPtHF",  "fake",  "miss",  "JER",  "PUprof"};





// Function prototype for flexible-scale function 
double Function_Mu(double s1, double s2 );

double Function_q(double q, double pt );
double Function_pt(double q, double pt );


TString rn() {return Form("%d",rand());}

using namespace std;


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
        if(!infile.good()) {
            cout << "File " << fName <<" does not exist." << endl;
            exit(1);
        }
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
            p.errUnc  /= 100;

            //p.sigma *= 0.97; //RADEK test

            //if(abs(p.yMin - 0) < 0.1) //first y-bin 0.5% //RADEK change
                //p.errUnc /= 2;

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
            //p.errs[1] = 0.04; //Radek test
            //cout << "HelenkaKarel " << p.errs[1] << endl;

            if(p.sigma > 0)
                dataNow.push_back(p);
            //cout << endl << endl;;
            // process pair (a,b)
        }
        assert(dataNow.size() > 10);
        return dataNow;

    }


    //Read theory histogram pdfName, as and scale variation s (the NP/EW corrections are applied)
    //tag for example 16ak4 or 16ak7
    static vector<vector<TH1D*>> readHistos(TString pdfName, double as, int s, TString tag, TString order = "nll")
    {

        int asI = round(as*1000);
        //Size of the histogram depends on the as value and #pdf for given pdf set
        
        vector<vector<TH1D*>> vTh; //indexes - [pdfVar][y]
        vTh.resize(1);
        if(pdfName.Contains("CT14")  && asI == 118)
            vTh.resize(56+1);
        else if(pdfName.Contains("HERAPDF")  && asI == 118)
            vTh.resize(28+1);
        else if(pdfName.Contains("NNPDF31")  && asI == 118)
            vTh.resize(100+1);
        else if(pdfName.Contains("ABMP16_5")  && asI == 118)
            vTh.resize(30);
        else if(pdfName.Contains("MMHT2014nnlo68cl")  && asI == 118)
            vTh.resize(50+1);

        int idCnt = 0;
        for(int ipdf = 0; ipdf < vTh.size(); ++ipdf) {
            for(int y = 0; y < 5; ++y) {
                TH1D *hTmp = (TH1D*) fTh->Get(pdfName+ Form("_y%d_as0%d_scale%d_pdf%d",y,asI, s, ipdf) );
                cout << pdfName+ Form("_y%d_as0%d_scale%d_pdf%d",y,asI, s, ipdf)  << endl;
                assert(hTmp);
                TH1D *h = (TH1D*) hTmp->Clone(rn()) ;
                vTh[ipdf].push_back(h);
            }
        }

        for(int ipdf = 0; ipdf < vTh.size(); ++ipdf) {
            for(int y = 0; y < 5; ++y) {
                applyNPEW(vTh[ipdf][y],  y, tag);

                TString tagN = tag;
                if(tag.Contains("ak4")) tagN = "_ak4";
                else if(tag.Contains("ak7")) tagN = "_ak7";
                else assert(0);

                if(order.Contains("nll")) applyKfactor(vTh[ipdf][y], y, "kFactorNLL"+tagN);
                else if(order.Contains("nnlo")) applyKfactor(vTh[ipdf][y], y, "kFactorNNLO"+tagN);
                //else if(order.Contains("nlo"))
                //applyKfactor(vTh[ipdf][y], y, "kFactorNNLO_ak4");
            }
        }

        return vTh;
    }

    //Read theory histograms for PDF pdfName (all alphaS (as) and all scale choices (s))
    void readAllTheory(TString pdfName, TString tag, TString order) {
        for(auto as: pdfAsVals.at(pdfName)) {
            cout << pdfName <<" "<< as << endl;
            int asI = round(as*1000);
            thHists[pdfName][asI].resize(7);
            for(int s = 0; s < 7; ++s)
                thHists[pdfName][asI][s] = readHistos(pdfName, as, s, tag, order);
        }
    }

    //Read theory histograms for PDF pdfName 
    void readSingleTheory(TString pdfName, TString tag, TString order) {
        for(auto as: pdfAsVals.at(pdfName)) {
            int asI = round(as*1000);
            if(asI != 118) continue;
            cout << pdfName <<" "<< as << endl;
            thHists[pdfName][asI].resize(1);
            for(int s = 0; s < 1; ++s)
                thHists[pdfName][asI][s] = readHistos(pdfName, as, s, tag, order);
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

            /*
            for(int i = 0; i < thHist118.size(); ++i) { //from 0?
                double thNom = thHist118[0][y]->GetBinContent(binId);
                double diff = (thHist118[i][y]->GetBinContent(binId) - thNom) / thNom;
                p.thErrs.push_back(diff);
            }
            */

            //Symetric hessian
            if(pdfName.Contains("ABMP16") || pdfName.Contains("NNPDF31")) {
                for(int i = 1; i < thHist118.size(); ++i) { 
                    double thNom = thHist118[0][y]->GetBinContent(binId);
                    double diff = (thHist118[i][y]->GetBinContent(binId) - thNom) / thNom;
                    p.thErrs.push_back(diff);
                }
            }
            //Assymetrick hessian - HERAPDF or CT14
            else {
                for(int i = 0; i < thHist118.size()/2; ++i) { 
                    double thNom = thHist118[0][y]->GetBinContent(binId);
                    double diff = (thHist118[2*i+1][y]->GetBinContent(binId) - thHist118[2*i+2][y]->GetBinContent(binId)) / thNom;

                    if(pdfName.Contains("CT14")) //CT14
                        p.thErrs.push_back(diff/2 * 1./1.645);
                    else
                        p.thErrs.push_back(diff/2);
                }
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

    void printHighest(const TVectorD &s, int n) {
        //TString sh;
        map<double, int> shifts;
        for(int i = 0; i < s.GetNrows(); ++i)
            shifts[-abs(s(i))] = i;
        int i = 0;
        for(auto el : shifts) {
            if(i < n) {
                if(el.second < ErrNames.size())
                    cout << ErrNames[el.second] <<" : "<< s(el.second) << ", ";
                else
                    cout << el.second <<" : "<< s(el.second) << ", ";
            }
            ++i;
        }
        cout << endl;
    }


    double getChi2All()
    {
        TVectorD s = getShiftsAll();
        //s.Print();
        //printHighest(s, 4);
        //s.Reset();
        //for(int i = 0; i < ErrNames.size(); ++i) s[i] = 0;

        double chi2 = getChi2All(s);
        
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

            double ref = p.sigma;
            double th  = p.th;
            double C = pow(p.sigma*p.errStat,2) + pow(p.sigma*p.errUnc,2);
            for(int j = 0; j < p.errs.size(); ++j)
                for(int k = 0; k < p.errs.size(); ++k)
                    mat(j,k) += 1./C * ref*ref * p.errs[j]*p.errs[k];

            for(int j = 0; j < p.errs.size(); ++j)
                yVec(j) += - 1./C * (p.sigma - th) * ref * p.errs[j];

        }

        for(int j = 0; j < data[0].errs.size(); ++j)
            mat(j,j) += 1;

        //Solve 
        TDecompSVD svd(mat);
        Bool_t ok;
        const TVectorD sh = svd.Solve(yVec, ok);

        return sh;
    }

    //Get the vector with the nuissence parameters, including theor unc. (values which minimize chi2)
    TVectorD getShiftsAll()
    {
        int nErr = data[0].errs.size() + data[0].thErrs.size();
        TMatrixD mat(nErr, nErr);
        TVectorD yVec(nErr);
        //Calculate the optimal shifts
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double th  = p.th;
            double ref = p.sigma;
            double C = pow(p.sigma*p.errStat,2) + pow(p.sigma*p.errUnc,2);
            for(int j = 0; j < nErr; ++j)
            for(int k = 0; k < nErr; ++k) {
                double thJ = (j < p.errs.size()) ? p.errs[j] : p.thErrs[j-p.errs.size()];
                double thK = (k < p.errs.size()) ? p.errs[k] : p.thErrs[k-p.errs.size()];
                mat(j,k) += 1./C * ref*ref * thJ*thK;
            }

            for(int j = 0; j < nErr; ++j) {
                double thJ = (j < p.errs.size()) ? p.errs[j] : p.thErrs[j-p.errs.size()];
                yVec(j) += - 1./C * (p.sigma - th) * ref * thJ;
            }

        }

        for(int j = 0; j < nErr; ++j)
            mat(j,j) += 1;

        //Solve 
        TDecompSVD svd(mat);
        Bool_t ok;
        const TVectorD sh = svd.Solve(yVec, ok);

        return sh;
    }

    // HERA chi2 fit with theory unc
    // http://www-h1.desy.de/psfiles/papers/desy15-039.pdf
    TVectorD getShiftsHERAall()
    {
        int nErr = data[0].errs.size() + data[0].thErrs.size();
        TMatrixD mat(nErr, nErr);
        TVectorD yVec(nErr);
        //Calculate the optimal shifts
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double m  = p.th;
            double mu = p.sigma;
            double C =  m*mu*pow(p.errStat,2) + m*m*pow(p.errUnc,2);
            for(int j = 0; j < nErr; ++j)
            for(int k = 0; k < nErr; ++k) {
                double thJ = (j < p.errs.size()) ? p.errs[j] : p.thErrs[j-p.errs.size()];
                double thK = (k < p.errs.size()) ? p.errs[k] : p.thErrs[k-p.errs.size()];
                mat(j,k) += 1./C * m*m * thJ*thK;
            }

            for(int j = 0; j < nErr; ++j) {
                double thJ = (j < p.errs.size()) ? p.errs[j] : p.thErrs[j-p.errs.size()];
                yVec(j) += - 1./C * (p.sigma - m) * m * thJ;
            }

        }

        for(int j = 0; j < nErr; ++j)
            mat(j,j) += 1;

        //Solve 
        TDecompSVD svd(mat);
        Bool_t ok;
        const TVectorD sh = svd.Solve(yVec, ok);

        return sh;
    }


    // HERA chi2 fit with theory unc (with fixed shift)
    // http://www-h1.desy.de/psfiles/papers/desy15-039.pdf
    // iShift - idOf the fixed shift, shVal - its val
    TVectorD getShiftsHERAall(int iShift, double shVal)
    {
        int nErr = data[0].errs.size() + data[0].thErrs.size();
        int nErrN= nErr-1; //new number of entries

        //map: newIndex -> oldIndex
        vector<double> indxMap;
        for(int i = 0; i < nErr; ++i) {
            if(i == iShift) continue;
            indxMap.push_back(i);
        }
        
        TMatrixD mat(nErrN, nErrN);
        TVectorD yVec(nErrN);
        //Calculate the optimal shifts
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double m  = p.th;
            double mu = p.sigma;
            double C =  m*mu*pow(p.errStat,2) + m*m*pow(p.errUnc,2);
            for(int j = 0; j < nErrN; ++j)
            for(int k = 0; k < nErrN; ++k) {
                int jG = indxMap[j]; //to old index
                int kG = indxMap[k];
                double thJ = (jG < p.errs.size()) ? p.errs[jG] : p.thErrs[jG-p.errs.size()];
                double thK = (kG < p.errs.size()) ? p.errs[kG] : p.thErrs[kG-p.errs.size()];
                mat(j,k) += 1./C * m*m * thJ*thK;
            }

            for(int j = 0; j < nErrN; ++j) {
                int jG = indxMap[j];
                double thJ = (jG < p.errs.size()) ? p.errs[jG] : p.thErrs[jG-p.errs.size()];
                double thI = (iShift < p.errs.size()) ? p.errs[iShift] : p.thErrs[iShift-p.errs.size()];
                yVec(j) += - 1./C * (mu - m + shVal*thI*m) * m * thJ; //including the fixed shift
            }

        }

        for(int j = 0; j < nErrN; ++j)
            mat(j,j) += 1;

        //Solve 
        TDecompSVD svd(mat);
        Bool_t ok;
        const TVectorD sh = svd.Solve(yVec, ok);

        //Inser the fixed value to the shifts
        TVectorD shNew(nErr); 
        for(int i = 0; i < nErrN; ++i) {
            shNew(indxMap[i]) = sh(i);
        }
        shNew(iShift) = shVal;

        return shNew;
    }











    //get the chi2 vale, the nuisence vector is as an input
    double getChi2(const TVectorD &s)
    {
        assert(data[0].errs.size() == s.GetNrows());

        //Evaluate the chi2
        double chi2 = 0;
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double th = p.th;
            double corErr = 0;
            for(int i = 0; i < p.errs.size(); ++i)
                corErr += s(i) * p.errs[i];

            double C = pow(p.sigma*p.errStat,2) + pow(p.sigma*p.errUnc,2);
            chi2 += pow(p.sigma - th + p.sigma*corErr, 2) / C;
        }

        for(int j = 0; j < data[0].errs.size(); ++j)
            chi2 += pow(s(j),2);

        return chi2;
    }

    //get the chi2 value, the nuisence vector is as an input, theor unc included
    double getChi2All(const TVectorD &s)
    {
        assert(data[0].errs.size() + data[0].thErrs.size() == s.GetNrows());

        //Evaluate the chi2
        double chi2 = 0;
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double th  = p.th;
            double ref = p.sigma;
            double corErr = 0;
            for(int i = 0; i < p.errs.size() + p.thErrs.size(); ++i) {
                double thI = (i < p.errs.size()) ? p.errs[i] : p.thErrs[i-p.errs.size()];
                corErr += s(i) * thI;
            }

            double C = pow(p.sigma*p.errStat,2) + pow(p.sigma*p.errUnc,2);
            chi2 += pow(p.sigma - th  + corErr*ref, 2) / C;
        }

        for(int j = 0; j < data[0].errs.size() + data[0].thErrs.size(); ++j)
            chi2 += pow(s(j),2);

        return chi2;
    }


    //get the chi2 value, the nuisence vector is as an input, theor unc included
    //Hera furmula http://www-h1.desy.de/psfiles/papers/desy15-039.pdf
    double getChi2HERAall(const TVectorD &s)
    {
        assert(data[0].errs.size() + data[0].thErrs.size() == s.GetNrows());

        //Evaluate the chi2
        double chi2 = 0;
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double m   = p.th;
            double mu  = p.sigma;
            //double th  = p.th;
            //double ref = p.sigma;
            double corErr = 0;
            for(int i = 0; i < p.errs.size() + p.thErrs.size(); ++i) {
                double thI = (i < p.errs.size()) ? p.errs[i] : p.thErrs[i-p.errs.size()];
                corErr += s(i) * thI;
            }

            double C = m*mu*pow(p.errStat,2) + m*m*pow(p.errUnc,2);
            //chi2 += pow(mu - m  + corErr*m, 2) / C;
            chi2 += pow(m   - corErr*m  - mu, 2) / C;

            //Log penalty
            chi2 += log(C / ((pow(p.errStat,2)+pow(p.errUnc,2))*mu*mu));
        }

        for(int j = 0; j < data[0].errs.size() + data[0].thErrs.size(); ++j)
            chi2 += pow(s(j),2);

        return chi2;
    }


    //With theory, but without correlated part
    pair<double,double> getChi2HERAallPartial(const TVectorD &s)
    {
        assert(data[0].errs.size() + data[0].thErrs.size() == s.GetNrows());

        //Evaluate the chi2
        double chi2Lin = 0;
        double chi2Log = 0;
        for(const auto &p : data) {
            if(!Cut(p)) continue;

            double m   = p.th;
            double mu  = p.sigma;
            //double th  = p.th;
            //double ref = p.sigma;
            double corErr = 0;
            for(int i = 0; i < p.errs.size() + p.thErrs.size(); ++i) {
                double thI = (i < p.errs.size()) ? p.errs[i] : p.thErrs[i-p.errs.size()];
                corErr += s(i) * thI;
            }

            double C = m*mu*pow(p.errStat,2) + m*m*pow(p.errUnc,2);
            //chi2 += pow(mu - m  + corErr*m, 2) / C;
            chi2Lin += pow(m   - corErr*m  - mu, 2) / C;

            //Log penalty
            chi2Log += log(C / ((pow(p.errStat,2)+pow(p.errUnc,2))*mu*mu));
        }

        //for(int j = 0; j < data[0].errs.size() + data[0].thErrs.size(); ++j)
            //chi2 += pow(s(j),2);

        return {chi2Lin, chi2Log};
    }







    //Get chi2 based on covariance matrix
    double getChi2cov(vector<int> indx)
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
        //for(int k = 0; k < dataF[0].errs.size(); ++k) {
        for(int k : indx) {
            for(int i = 0; i < dataF.size(); ++i) 
            for(int j = 0; j < dataF.size(); ++j) {
                //if(dataF[i].yMin != dataF[j].yMin)// && k <= 15) continue;
                //if(i != j) continue;
                CovSys(i,j) += dataF[i].errs[k]*dataF[j].errs[k]  *  dataF[i].sigma * dataF[j].sigma;
            }
        }

        //Fill pdf sys matrix

        double Ccorr = 1;
        /*
        if(dataF[0].thErrs.size() == 28) //CT14
            Ccorr = 1*1./(1.64*1.64);
        else if(dataF[0].thErrs.size() == 29) //HERAPDF
            Ccorr = 1*1./(2);
        else if(dataF[0].thErrs.size() == 101) //NNPDF31
            Ccorr = 1*1./(2);
        else {
            cout << "pdf err size " << dataF[0].thErrs.size()  << endl;
            assert(0);
        }
        */

        TMatrixD CovPDF(dataF.size(), dataF.size());
        for(int k = 0; k < dataF[0].thErrs.size(); ++k) {
            for(int i = 0; i < dataF.size(); ++i) 
            for(int j = 0; j < dataF.size(); ++j) {
                //if(dataF[i].yMin != dataF[j].yMin) continue;
                CovPDF(i,j) += dataF[i].thErrs[k]*dataF[j].thErrs[k]  *  dataF[i].sigma * dataF[j].sigma * Ccorr;
            }
        }

        //TMatrixD Cov = CovStat + CovSys;// +  CovPDF;// + CovSys + CovPDF;
        TMatrixD Cov = CovStat + CovSys +  CovPDF;// + CovSys + CovPDF;

        
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

        double chi2 = 0;
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

        double chi2N = getChi2All();

        auto shifts  = getShiftsHERAall();
        double chi2H = getChi2HERAall(shifts);


        vector<int> indx;
        for(int i = 0; i < data[0].errs.size(); ++i) {
            if(i != -1) indx.push_back(i); //remove luminosity
        }
        double chi2C = getChi2cov(indx);
        cout <<"chi2 "<< chi2N <<" "<< chi2C << " "<< chi2H << endl; 
        return chi2C;
        //return getChi2();
    }


    void ScanChi2(TString pdfName, int scale = 0) {

        for(auto as: pdfAsVals.at(pdfName)) {
            fillTheory(pdfName, as, scale);

            for(int i = 0; i < data[0].errs.size(); ++i) {
                vector<int> indx;
                indx.push_back(i);
                double chi2 =  getChi2cov(indx);
                cout << as <<" "<< i <<" "<< chi2 << endl;
            }
        }
        //return getChi2();
    }




    //print the chi2 table
    void printChi2Table()
    {
        vector<TString> pdfSets = {"CT14nnlo", "HERAPDF20_NNLO", "NNPDF31_nnlo", "ABMP16_5_nnlo"};

        for(auto p : pdfSets)
            cout << p << " & ";
        cout << endl;

        for(int y = 0; y < 4; ++y) { //over rapidities
            cout << y*0.5 <<" & " << (y+1)*0.5 << " & ";
            Cut = [y](point p) { return ( abs(y*0.5-p.yMin) < 0.1 &&  p.sigma != 0);};
            int ndf = getNpoints();

            cout << ndf << " & ";
            for(auto pdfSet : pdfSets) {
                double chi2now = calcChi2(pdfSet, 0.118);
                //cout  <<chi2now << " / " << ndf << " " ;
                cout  <<chi2now << " & ";
            }
            cout << "//" << endl;

        }

        
        { //Total chi2 
            cout << "Total &&";
            Cut = [](point p) { return ( abs(p.yMin) < 1.6 &&  p.sigma != 0);};
            int ndf = getNpoints();

            cout << ndf << " & ";
            for(auto pdfSet : pdfSets) {
                double chi2now = calcChi2(pdfSet, 0.118);
                cout  <<chi2now <<  " & " ;
            }
            cout << "//" <<  endl;

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

			
    //
    void TheoryPlotter()
    {	



    }

    void plotReview(TString pdfName, double as, int scale = 0)
    {
        const int nYbins = 4;
        Cut = [](point p) { return ( abs(p.yMin) < 0.5*(nYbins-0.9) &&  p.sigma != 0 && p.ptMin > 95);};

        gStyle->SetOptStat(0);
        fillTheory(pdfName, as, scale);


        int nSys = data[0].errs.size();
        int nTh  = data[0].thErrs.size();

        auto shifts  = getShiftsHERAall();
        double chi2H = getChi2HERAall(shifts);

        //Uncertainties of the shifts
        TVectorD shiftsUnc(shifts.GetNrows());
        for(int i = 0; i < shifts.GetNrows(); ++i) {
            auto shiftsU  = getShiftsHERAall(i, shifts(i)+1); //up-variation
            double dChi2 = getChi2HERAall(shiftsU) - chi2H;
            shiftsUnc(i) = 1./sqrt(dChi2);
        }
        

        /*
        auto shiftsT   = getShiftsHERAall(0, shifts(0));
        auto shiftsTU  = getShiftsHERAall(0, shifts(0)+1);
        auto shiftsTD  = getShiftsHERAall(0, shifts(0)-1);
        for(int i = 0; i < shifts.GetNrows(); ++i) {
            cout << i <<" "<< shifts(i) <<" "<<shiftsT(i) << endl;
        }
        cout << "Chi2 comparison " << getChi2HERAall(shifts) <<" : "<< getChi2HERAall(shiftsT) << " "<< getChi2HERAall(shiftsTU)<<" "<< getChi2HERAall(shiftsTD) << endl;
        exit(0);
        */

        int ndfT  = getNpoints();
        assert(shifts.GetNrows() == nSys + nTh);


        //Get chi2 for all rap bins

        vector<double> chi2NY(nYbins), chi2LY(nYbins);
        vector<int> ndfY(nYbins);
        for(int y = 0; y < nYbins; ++y) {
            Cut = [y](point p) { return ( round(abs(2*p.yMin)) == y &&  p.sigma != 0 && p.ptMin > 95);};
            tie(chi2NY[y],chi2LY[y]) = getChi2HERAallPartial(shifts);
            ndfY[y]  = getNpoints();
        }





        //Shifts to histogram
        TH1D *hShifts = new TH1D(rn(), "", shifts.GetNrows(), 0.5, shifts.GetNrows()+0.5);
        for(int i = 0; i < shifts.GetNrows(); ++i) {
            hShifts->SetBinContent(i+1, shifts[i]);
            hShifts->SetBinError(i+1, shiftsUnc[i]);
            if(i < nSys) hShifts->GetXaxis()->SetBinLabel(i+1, ErrNames[i]);
        }




        //Get binning of the data
        vector<vector<double>> bins(nYbins);
        for(auto p : data) {
            int y = round(p.yMin*2);
            if(y >= nYbins || p.ptMin < 95) continue;
            bins[y].push_back(p.ptMin);
            bins[y].push_back(p.ptMax);
        }
        for(auto &bin : bins) {
            sort(bin.begin(), bin.end());
            auto it = unique(bin.begin(), bin.end());
            bin.resize(distance(bin.begin(), it));
        }

        //Init all histograms
        vector<TH1D*> hData(nYbins);
        vector<TH1D*> hTh(nYbins);
        vector<TH1D*> hThShTot(nYbins);
        vector<vector<TH1D*>> hShData(nSys);
        vector<vector<TH1D*>> hShTh(nTh);

        for(int y = 0; y < nYbins; ++y) {
           hData[y] = new TH1D(rn(), "", bins[y].size()-1, bins[y].data()); 
           hTh[y]   = new TH1D(rn(), "", bins[y].size()-1, bins[y].data()); 
           hThShTot[y] = new TH1D(rn(), "", bins[y].size()-1, bins[y].data()); 
           for(int s = 0; s < nSys; ++s) {
               hShData[s].resize(nYbins);
               hShData[s][y] = new TH1D(rn(), "", bins[y].size()-1, bins[y].data()); 
           }

           for(int s = 0; s < nTh; ++s) {
               hShTh[s].resize(nYbins);
               hShTh[s][y] = new TH1D(rn(), "", bins[y].size()-1, bins[y].data()); 
           }
        }





        //Fill theory and data
        for(auto p : data) {
            int y = round(p.yMin*2);
            if(y >= nYbins) continue;
            int ipt = hData[y]->FindBin((p.ptMin+p.ptMax)/2);
            hData[y]->SetBinContent(ipt, p.sigma);
            hData[y]->SetBinError(ipt, hypot(p.errStat,p.errUnc)*p.sigma);
            hTh[y]->SetBinContent(ipt, p.th);
            hTh[y]->SetBinError(ipt, 0);


            //Fill shifted theory
            double shTot = 0;
            for(int s = 0; s < p.errs.size(); ++s) {
                shTot += p.errs[s]*shifts[s]*p.th;
                hShData[s][y]->SetBinContent(ipt, -p.errs[s]*shifts[s]);
                hShData[s][y]->SetBinError(ipt, 0);
            }
                
            for(int s = 0; s < p.thErrs.size(); ++s) {
                shTot += p.thErrs[s]*shifts[p.errs.size()+s]*p.th;
                hShTh[s][y]->SetBinContent(ipt, -p.thErrs[s]*shifts[p.errs.size()+s]);
                hShTh[s][y]->SetBinError(ipt, 0);
            }
            hThShTot[y]->SetBinContent(ipt, p.th-shTot);
            hThShTot[y]->SetBinError(ipt, 0);

        }
        

        //Get important data shifts
        map<double,int> shImportance;
        for(int s = 0; s < nSys; ++s) {
            double M = -1e10;
            for(int y = 0; y < nYbins; ++y) {
                M = max(M, abs(hShData[s][y]->GetBinContent(2)));
                M = max(M, abs(hShData[s][y]->GetBinContent(5)));
                M = max(M, abs(hShData[s][y]->GetBinContent(10)));
            }
            shImportance[M] = s;
        }
        auto it = shImportance.begin();
        cout << "Radek start " << shImportance.size() <<" "<< nSys <<endl;
        std::advance(it, shImportance.size()-5);
        shImportance.erase(shImportance.begin(), it);
        cout << "Radek after " << shImportance.size() << endl;
        map<int,int> shImp;
        vector<int> myCols = {kBlue, kGreen+2, kYellow+2, kViolet, kCyan+2,
                         kPink+6, kOrange+3, kAzure-4, kGray+3, kGreen-6};
        int ic = 0;
        for(auto s : shImportance) {
            //cout << "Radek insering " << s.second << endl;
            shImp[s.second] = myCols[ic++];
        }



        TCanvas *can = new TCanvas(rn(), "", 1000, 700);
        SetTopBottom(0.05, 0.15);
        DividePad({1}, {1,1,1,1,1});
        //DivideTransparent({1}, {1,0,1,0,1,0,1,0.2,1});

        //Ratio to theory
        can->cd(1);
        DividePad(vector<double>(nYbins,1.), {1});
        for(int y = 0; y < nYbins; ++y) {
            can->cd(1)->cd(y+1);
            gPad->SetLogx();
            auto hDataR = (TH1D*) hData[y]->Clone(rn());
            auto hThShR = (TH1D*) hThShTot[y]->Clone(rn());
            auto hThR   = (TH1D*) hTh[y]->Clone(rn());

            //int ipt = hDataR->FindBin(240);
            //cout << "Data: " << hDataR->GetBinContent(ipt) << endl;
            //cout << "NNLO: " << hThR->GetBinContent(ipt) << endl;
            //exit(0);

            hDataR->Divide(hTh[y]);
            hThShR->Divide(hTh[y]);
            hThR->Divide(hTh[y]);

            hDataR->SetLineColor(kBlack);
            hThShR->SetLineColor(kBlue);
            hThR->SetLineColor(kRed);

            hThR->Draw();
            hThShR->Draw("same ][");
            hDataR->Draw("same");
            GetYaxis()->SetRangeUser(0.6,1.5);
            GetYaxis()->SetNdivisions(404);
            GetYaxis()->SetTitle("data/theory");
            SetFTO({20}, {10}, {1.3, 1.5, 0.4, 3.4});

            DrawLatexUp(-1.05, yBins[y]);
        }

        can->cd(1)->cd(4);
        DrawLatexUp(1.3, Form("#alpha_{S} = %.3f : #chi^{2} = %.1f / %d", as, chi2H, ndfT), 20, "c");


        //Ratio to shifted theory
        can->cd(2);
        DividePad(vector<double>(nYbins,1.), {1});
        for(int y = 0; y < nYbins; ++y) {
            can->cd(2)->cd(y+1);
            gPad->SetLogx();
            auto hDataR = (TH1D*) hData[y]->Clone(rn());
            auto hThShR = (TH1D*) hThShTot[y]->Clone(rn());
            hDataR->Divide(hThShTot[y]);
            hThShR->Divide(hThShTot[y]);

            hDataR->SetLineColor(kBlack);
            hThShR->SetLineColor(kBlue);
            hThShR->Draw("][");
            hDataR->Draw("same");
            GetYaxis()->SetRangeUser(0.9,1.1);
            GetYaxis()->SetTitle("data/theory'");
            GetYaxis()->SetNdivisions(404);

            DrawLatexUp(-1.2, Form("#chi^{2} = %.1f / %d", chi2NY[y]+ chi2LY[y], ndfY[y]),20);
            SetFTO({20}, {10}, {1.3, 1.5, 0.4, 3.4});


        }

        //Shifts pt-spectra review
        can->cd(3);
        DividePad(vector<double>(nYbins,1.), {1});
        for(int y = 0; y < nYbins; ++y) {
            can->cd(3)->cd(y+1);
            gPad->SetLogx();


            auto hThShR = (TH1D*) hThShTot[y]->Clone(rn());
            hThShR->Divide(hTh[y]);
            for(int k = 1; k < hThShR->GetNbinsX(); ++k)
                hThShR->SetBinContent(k, hThShR->GetBinContent(k)-1);

            hThShR->SetLineColor(kBlack);
            hThShR->Draw("hist");


            for(int s = 0; s < nTh; ++s) {
                hShTh[s][y]->SetLineColor(kRed);
                hShTh[s][y]->Draw("same ][");
            }

            for(int s = 0; s < nSys; ++s) {
                if(shImp.count(s)) continue;
                hShData[s][y]->SetLineColor(kOrange);
                hShData[s][y]->Draw("same ][");
            }
            for(auto sEl : shImp) {
                int s = sEl.first;
                int c = sEl.second;
                hShData[s][y]->SetLineColor(c);
                hShData[s][y]->Draw("same ][");
            }
            hThShR->Draw("same");


            GetYaxis()->SetRangeUser(-0.2,0.2);
            GetYaxis()->SetNdivisions(404);
            GetYaxis()->SetTitle("Rel. Unc.");
            SetFTO({20}, {10}, {1.3, 1.5, 0.4, 3.4});


            auto leg = newLegend(kPos1);
            leg->AddEntry(hThShR, "total shift", "l");
            DrawLegends({leg}, true);
        }


        //MAIN shifts pt-spectra review
        can->cd(4);
        DividePad(vector<double>(nYbins,1.), {1});
        for(int y = 0; y < nYbins; ++y) {
            can->cd(4)->cd(y+1);
            gPad->SetLogx();

            //Total shift
            auto hThShR = (TH1D*) hThShTot[y]->Clone(rn());
            hThShR->Divide(hTh[y]);
            for(int k = 1; k < hThShR->GetNbinsX(); ++k)
                hThShR->SetBinContent(k, hThShR->GetBinContent(k)-1);


            //Total PDF shift
            auto hPDFTotSh = (TH1D*)hShTh[0][y]->Clone();
            for(int s = 1; s < nTh; ++s) {
                hPDFTotSh->Add(hShTh[s][y]);
            }

            //Total Sys shift (exclusing dominant)
            auto hDataTotSh = (TH1D*)hShData[0][y]->Clone();
            hDataTotSh->Reset();
            for(int s = 0; s < nSys; ++s) {
                if(shImp.count(s)) continue;
                hDataTotSh->Add(hShData[s][y]);
            }


            hThShR->Draw("hist");



            //Main shifts
            for(auto sEl : shImp) {
                int s = sEl.first;
                int c = sEl.second;
                hShData[s][y]->SetLineColor(c);
                hShData[s][y]->Draw("same ][");
            }

            hPDFTotSh->SetLineColor(kRed);
            hPDFTotSh->Draw("same ][");

            hDataTotSh->SetLineColor(kOrange);
            hDataTotSh->Draw("same ][");


            hThShR->SetLineColor(kBlack);
            hThShR->Draw("same ][");

            GetYaxis()->SetRangeUser(-0.2,0.2);
            GetYaxis()->SetNdivisions(404);
            GetYaxis()->SetTitle("Rel. Unc.");
            SetFTO({20}, {10}, {1.3, 1.5, 0.4, 3.4});
        }


        //Shifts
        can->cd(5);

        double chiSys = 0, chiTh = 0;
        for(int i = 0; i < nSys; ++i)
            chiSys += pow(shifts[i],2);
        for(int i = 0; i < nTh; ++i)
            chiTh  += pow(shifts[nSys+i],2);

        hShifts->SetLineColor(kBlack);
        hShifts->Draw();

        //Draw important shifts with corresponding collors
        for(auto s : shImp) {
            auto hNow = (TH1D*) hShifts->Clone(rn());
            for(int i = 1; i <= hNow->GetNbinsX(); ++i) {
                if(i-1 != s.first) hNow->SetBinContent(i, -1000);
            }
            hNow->SetLineColor(s.second);
            hNow->Draw("same");
        }

        SetFTO({20}, {10}, {1.3, 1.5, 0.4, 2.8});
        GetYaxis()->SetRangeUser(-2.9, 3.6);
        GetYaxis()->SetTitle("shift");

        TLine *l = new TLine;
        l->SetLineStyle(2);
        l->DrawLine(nSys+0.5, -3, nSys+0.5, 4);
        l->DrawLine(1+0.5, -3, 1+0.5, 4);
        l->DrawLine(0.5, -1, nSys+nTh+0.5, -1);
        l->DrawLine(0.5,  1, nSys+nTh+0.5,  1);
        DrawLatexUp(-1, Form("    Data: #chi^{2} = %.1f / %d", chiSys, nSys), -1, "l");
        DrawLatexUp(-1, Form("PDF: #chi^{2} = %.1f / %d   ", chiTh, nTh), -1, "r");
        GetYaxis()->SetNdivisions(404);



        int asI = round(1000*as);
        int asF = round(1000* pdfAsVals.at(pdfName).front());
        int asB = round(1000* pdfAsVals.at(pdfName).back());
        if(asI == asF)
            can->SaveAs("rew.pdf(");
        else if(asI == asB)
            can->SaveAs("rew.pdf)");
        else
            can->SaveAs("rew.pdf");
    }

};













int main(int argc, char** argv)
{
    fTh  = TFile::Open("cmsJetsAsScan_ak7.root");  //NLO predictions

	asFitter asfit;
    //asfit.data = asfit.readData("xFitterTables/patrick16ak4.txt");
    asfit.data = asfit.readData("xFitterTables/data16ak7NewNew.txt");
    //asfit.readAllTheory("CT14nnlo");
    //asfit.readSingleTheory("HERAPDF20_NNLO");

    //TString curPDF = "NNPDF31_nnlo";
    TString curPDF = "CT14nnlo";
    //TString curPDF = "ABMP16_5_nnlo";
    //TString curPDF = "HERAPDF20_NNLO";

    asfit.readAllTheory(curPDF, "16ak7", "nnlo");
    //asfit.readSingleTheory("ABMP16_5_nnlo");
    //asfit.readSingleTheory("MMHT2014nnlo68cl");

    //asfit.printChi2Table();

    //return 0;
    //asfit.getAllChi2s();
    //return 0;

    //asfit.Cut = [](point p) { return ( abs(p.yMin) < 1.7 &&  p.sigma != 0 && p.ptMin > 96);};
    //asfit.ScanChi2("CT14nnlo");
    //return 0;


    cout << "Reading finished " << endl;
    asfit.Cut = [](point p) { return ( abs(p.yMin) < 1.7 &&  p.sigma != 0 && p.ptMin > 96);};


    for(auto as: pdfAsVals.at(curPDF))
        asfit.plotReview(curPDF, as, 0);
    //asfit.plotReview(curPDF, 0.118, 0);
    return 0;

    for(auto as: pdfAsVals.at(curPDF)) {
        //if(round(1000*as) != 118) continue;
        double chi2 = asfit.calcChi2(curPDF, as, 0);
        cout << "Helenka " << as <<" "<< chi2  <<" : "<< asfit.getNpoints() << endl;
    }

    //return 0;


    for(int y = 0; y < 4; ++y) {
        asfit.Cut = [y](point p) { return ( abs(y*0.5-p.yMin) < 0.1 &&  p.sigma != 0 && p.ptMin > 96);};
        for(auto as: pdfAsVals.at(curPDF)) {
            //if(round(1000*as) != 118) continue;
            double chi2 = asfit.calcChi2(curPDF, as, 0);
            cout << "Helenka " << y <<" "<<as <<" "<<  chi2  <<" : "<< asfit.getNpoints() << endl;
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
