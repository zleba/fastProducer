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

double Function_q(double q, double pt );
double Function_pt(double q, double pt );


TString rn() {return Form("%d",rand());}

using namespace std;

//All functions to have a list
vector<TH1D*> readHisto(fastNLOAlphas &fnlo);
vector<vector<TH1D*>> getScaleuncHistos(fastNLOAlphas &fnlo);
vector<vector<TH1D*>> getAsScaleuncHistos(TString pdfName, int R);
vector<vector<TH1D*>> getPDFuncHistos(fastNLOAlphas &fnlo);
//vector<vector<TH1D*>> getAsHistos(fastNLOAlphas &fnlo);
TH1D *rebin(TH1D *h, TH1D *hTemp);
void printHisto(TH1D *h);
void SaveHistos(vector<vector<TH1D*>> hist,  TString tag);
void SaveHistosByTitle(vector<vector<TH1D*>> hist);
void scanAsToFile(int R);
vector<vector<vector<TH1D*>>> calcXsections(int R, TString pdfName);

TString getLHAname(TString pdfName, int asI);



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

    vector<vector<double>> scales = { { 1, 1},
                              { 2, 2},
                              { 0.5, 0.5},
                              { 1, 2},
                              { 1, 0.5},
                              { 2, 1},
                              { 0.5, 1} };


    vector<vector<TH1D*>> histos;
    for(auto  s : scales) {
        fnlo.SetScaleFactorsMuRMuF(s[0], s[1]);
        auto hh = readHisto(fnlo);
        //cout << "RAdek before " << hh.size() << endl;
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


TString getLHAname(TString pdfName, int asI)
{
    TString tag1, tag2;
    if(pdfName == "CT14nlo") {
        tag1 = "CT14nlo_as_0";
    } else if(pdfName == "CT14nnlo") {
        tag1 = "CT14nnlo_as_0";

    } else if(pdfName == "HERAPDF20_NLO") {
        tag1 = "HERAPDF20_NLO_ALPHAS_";
    } else if(pdfName == "HERAPDF20_NNLO") {
        tag1 = "HERAPDF20_NNLO_ALPHAS_";

    } else if(pdfName == "NNPDF31_nlo") {
        tag1 = "NNPDF31_nlo_as_0";
    } else if(pdfName == "NNPDF31_nnlo") {
        tag1 = "NNPDF31_nnlo_as_0";

    } else if(pdfName == "ABMP16_5_nlo") {
        tag1 = "ABMP16als";
        tag2 = "_5_nlo";
    } else if(pdfName == "ABMP16_5_nnlo") {
        tag1 = "ABMP16als";
        tag2 = "_5_nnlo";
    } else if(pdfName == "MMHT2014nlo68cl" || "MMHT2014nnlo68cl") {
        return pdfName;
    }
    else
        exit(1);

    TString whole = tag1 + Form("%d", asI) + tag2;

    if(asI == 118) {
        if(pdfName.Contains("CT14")) 
            whole = pdfName;
        else if(pdfName.Contains("HERAPDF20_"))
            whole = pdfName + "_EIG";
        else if(pdfName.Contains("NNPDF31_nlo"))
            whole = "NNPDF31_nlo_as_0118_hessian";
        else if(pdfName.Contains("NNPDF31_nnlo"))
            whole = "NNPDF31_nnlo_as_0118_hessian";
    }
    if(pdfName.Contains("MMHT2014"))
        whole = pdfName;

    return whole;
}




//R = 4 or R = 7
vector<vector<TH1D*>> getAsScaleuncHistos(TString pdfName, int R)
{
    //fnlo.SetLHAPDFMember(0);

    vector<vector<double>> scales = { { 1, 1},
                              { 2, 2},
                              { 0.5, 0.5},
                              { 1, 2},
                              { 1, 0.5},
                              { 2, 1},
                              { 0.5, 1} };

    /*
    map<TString, vector<double> > pdfAsVals;
    pdfAsVals["CT14nlo"] = getRange(0.111, 0.123);
    pdfAsVals["CT14nnlo"] = pdfAsVals["CT14nlo"];

    pdfAsVals["HERAPDF20_NLO"] = getRange(0.111, 0.123);
    pdfAsVals["HERAPDF20_NNLO"] = pdfAsVals["HERAPDF20_NLO"];

    pdfAsVals["NNPDF31_nnlo"] = {0.112, 0.114, 0.116, 0.117, 0.118, 0.119, 0.120,  0.122};
    */


    cout << "Radek " << endl;
    vector<vector<TH1D*>> histos;
    for(double as : pdfAsVals.at(pdfName)) {
        int asI = round(as * 1000);

        TString whole = getLHAname(pdfName, asI);

        TString inFile;
        if(R == 4)      inFile = "theorFiles/InclusiveNJets_fnl5362h_v23_fix.tab";
        else if(R == 7) inFile = "theorFiles/InclusiveNJets_fnl5332h_v23_fix.tab";
        else assert(1);

        fastNLOAlphas  fnlo(inFile.Data(), whole.Data(), 0);
        //fastNLOAlphas  fnlo("theorFiles/suman/Fnlo_AK7_Eta1.tab", whole.Data(), 0);
        fnlo.SetAlphasMz(as, true);

        //cout << "Helenka " << as << endl;

        int sId = 0;
        for(auto  s : scales) {
            int nPDFs = fnlo.GetNPDFMembers();
            cout << sId << " "<< pdfName<<" : "<< whole <<" "<< nPDFs << endl;
            if(asI != 118)
                nPDFs = 1;
            for(int pdfId = 0; pdfId < nPDFs; ++pdfId) {
                fnlo.SetLHAPDFMember(pdfId);
                fnlo.SetScaleFactorsMuRMuF(s[0], s[1]);
                auto hh = readHisto(fnlo);

                for(int y = 0; y < hh.size(); ++y)
                    hh[y]->SetTitle(pdfName + Form("_y%d_as0%d_scale%d_pdf%d", y, asI, sId, pdfId));

                //cout << "RAdek before " << hh.size() << endl;
                histos.push_back(hh);
            }

            ++sId;
        }
        cout << "Helenka end" << endl;
    }
    //exit(0);

    return histos;
}




//Get histogram including up and dn pdf variation 
vector<vector<TH1D*>> getPDFuncHistos(fastNLOAlphas &fnlo)
{
    fnlo.SetLHAPDFMember(0);
    fnlo.SetScaleFactorsMuRMuF(1, 1);

    int nPDFs = fnlo.GetNPDFMembers();
    TString pdfName = fnlo.GetLHAPDFFilename();
    cout << pdfName << endl;
    double Fact = pdfName.Contains("CT14") ? 1.65 : 1;

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

//Get histo using interpolation
vector<TH1D*> getAsHisto(map<double, fastNLOAlphas*> &fnloMap, double as)
{
    fastNLOAlphas *fastL, *fastH;
    double asL, asH;
    for(auto fast : fnloMap) {
        if(fast.first <= as) {
            asL   = fast.first;
            fastL = fast.second;
        }
        if(fast.first >= as) {
            asH = fast.first;
            fastH = fast.second;
            break;
        }
    }
    auto hL = readHisto(*fastL);
    auto hH = readHisto(*fastH);

    vector<TH1D*> hAvg(hL.size());
    for(int y = 0; y < hL.size(); ++y) {
       hAvg[y] = (TH1D*) hL[y]->Clone(rn());
       hAvg[y]->Add(hL[y], hH[y],  (asH - as)/(asH-asL), (as - asL)/(asH-asL));
    }
    return hAvg;
}

//Get histogram including up and dn aS variation 
//input: fastNLO map ideally for 0.116, 0.117, 0.118, 0.119, 0.200
vector<vector<TH1D*>> getAsHistos(map<double, fastNLOAlphas*> &fnloMap)
{
    for(auto fast : fnloMap) {
        fast.second->SetLHAPDFMember(0);
        fast.second->SetScaleFactorsMuRMuF(1, 1);
        fast.second->SetAlphasMz(fast.first, true);
    }

    const double asErr = 0.0015;
    auto hU =  getAsHisto(fnloMap, 0.1180+asErr);
    auto hD =  getAsHisto(fnloMap, 0.1180-asErr);
    auto hC =  readHisto(*fnloMap.at(0.118));

    vector<TH1D*> hUp(hU.size());
    vector<TH1D*> hDn(hD.size());

    for(int y = 0; y < hC.size(); ++y) { //loop over y-bins

        hUp[y] = (TH1D*) hU[y]->Clone(rn());
        hDn[y] = (TH1D*) hD[y]->Clone(rn());

        for(int i = 1; i <= hC[y]->GetNbinsX(); ++i) { //loop over pt-bins
            double cnt = hC[y]->GetBinContent(i);
            double up = 0, dn = 0;
            up = max(cnt, max(hU[y]->GetBinContent(i),hD[y]->GetBinContent(i)));
            dn = min(cnt, min(hU[y]->GetBinContent(i),hD[y]->GetBinContent(i)));
            hUp[y]->SetBinContent(i, up);
            hDn[y]->SetBinContent(i, dn);
        }
    }
    return {hC, hUp, hDn};
}










//rebin differential histo according to template
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


void SaveHistosByTitle(vector<vector<TH1D*>> hist)
{
    for(int s = 0; s < hist.size(); ++s) {
        for(int y = 0; y < hist[s].size(); ++y) {
            //TString n = tag +"_"+ sysTag[s] +"_"+ Form("y%d", y);
            TString n = hist[s][y]->GetTitle();
            cout << "Saving to root " << s <<" "<< y <<" : "<<  n << endl;
            hist[s][y]->SetName(n);
            hist[s][y]->Write(n);
        }
    }
}


//Create the root file with many theoryes 
void scanAsToFile(int R)
{
    vector<TString> pdfList = {"CT14nlo", "CT14nnlo", "HERAPDF20_NLO", "HERAPDF20_NNLO",   "NNPDF31_nlo", "NNPDF31_nnlo", "ABMP16_5_nlo", "ABMP16_5_nnlo"};
    //vector<TString> pdfList = { "ABMP16_5_nlo", "ABMP16_5_nnlo"};


    TFile *fOut = new TFile(Form("cmsJetsAsScan_ak%d.root",R), "RECREATE");

    for(auto pdf : pdfList) {
        vector<vector<TH1D*>> histPDF    = getAsScaleuncHistos(pdf, R);
        SaveHistosByTitle(histPDF);
    }

    fOut->Write();
    fOut->Close();
}


void setNLO(fastNLOAlphas &fnlo)
{
    fnlo.SetContributionON(fastNLO::kFixedOrder,0,true);
    fnlo.SetContributionON(fastNLO::kFixedOrder,1,true);
    fnlo.SetUnits(fastNLO::kPublicationUnits);
}

vector<vector<vector<TH1D*>>> calcXsections(int R, TString pdfName)
{

    using namespace std;
    using namespace say;		// namespace for 'speaker.h'-verbosity levels
    using namespace fastNLO;	// namespace for fastNLO constants


    TString fastName;
    if(R == 4)      fastName = "theorFiles/InclusiveNJets_fnl5362h_v23_fix.tab";
    else if(R == 7) fastName = "theorFiles/InclusiveNJets_fnl5332h_v23_fix.tab";
    else assert(1);

    cout << "Helenka " << __LINE__ << endl;
    fastNLOAlphas fnlo(fastName.Data(), getLHAname(pdfName, 118).Data(), 0);
    cout << "Helenka " << __LINE__ << endl;
    setNLO(fnlo);
    cout << "Helenka " << __LINE__ << endl;

    vector<vector<TH1D*>> histScl = getScaleuncHistos(fnlo);
    cout << "Helenka " << __LINE__ << endl;
    vector<vector<TH1D*>> histPDF = getPDFuncHistos(fnlo);
    cout << "Helenka " << __LINE__ << endl;

    vector<vector<TH1D*>> histAs;
    if(!pdfName.Contains("MMHT2014")) {

        map<double,fastNLOAlphas*> fnloMap;
        //get alphaS unc vector
        vector<TString> names;

        for(auto as : pdfAsVals.at(pdfName)) {
        //for(int asI = 116; asI <= 120; ++asI) {
            int asI = round(as * 1000);
            if(asI < 116 || asI > 120) continue;
            TString pdfNameAs = getLHAname(pdfName, asI);
            cout << fastName.Data() <<" "<< pdfNameAs.Data() << endl;
            fnloMap.insert( make_pair(asI/1000., new fastNLOAlphas(fastName.Data(), pdfNameAs.Data(), 0)) );
            setNLO(*fnloMap.at(asI/1000.));
        }
        histAs = getAsHistos(fnloMap);
    }
    else { //dummy
        histAs = {histScl[0], histScl[0], histScl[0]}; 
    }

    return {histPDF, histScl, histAs};

    //fOut->Write();
    //fOut->Close();

}

void calcXsectionsAll(int R, vector<TString> pdfNames)
{
    vector<vector<vector<TH1D*>>>  histsPDF, histsScl, histsAs;

    for(auto pdf : pdfNames) {
        auto xsecs =  calcXsections(R, pdf);
        histsPDF.push_back(xsecs[0]);
        histsScl.push_back(xsecs[1]);
        histsAs. push_back(xsecs[2]);
    }

    TFile *fOut = new TFile(Form("theorFiles/cmsJetsNLO_AK%d.root", R), "RECREATE");
    for(int i = 0; i < histsPDF.size(); ++i) {
        SaveHistos(histsPDF[i], "hist"+  pdfNames[i] + "_PDF");
        SaveHistos(histsScl[i], "hist"+  pdfNames[i] + "_Scl");
        SaveHistos(histsAs[i],  "hist"+  pdfNames[i] + "_As");
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


    scanAsToFile(4);
    scanAsToFile(7);
    return 0;


    //calcXsectionsAll(4, {"CT14nlo",  "CT14nnlo", "HERAPDF20_NLO",  "HERAPDF20_NNLO", "NNPDF31_nlo", "NNPDF31_nnlo", "ABMP16_5_nlo", "ABMP16_5_nnlo", "MMHT2014nlo68cl", "MMHT2014nnlo68cl"});

    //TFile *fOut = new TFile("theorFiles/cmsJetsNLO_AK4.root", "RECREATE");
    ////calcXsections("theorFiles/fastnlo-cms-incjets-arxiv-1605.04436-xsec001.tab", "Old", "CT14nlo");
    //calcXsections(4, "CT14nlo");
    ////calcXsections("theorFiles/InclusiveNJets_fnl5362h_v23_fix.tab", "New", "CT14nlo");
    ////calcXsections(7, "CT14nlo");
    //fOut->Write();
    //fOut->Close();

    return 0;

    //vector<TH1D*> readHisto(fastNLOAlphas &fnlo, TString tabName, TString pdfName)

    //say::SetGlobalVerbosity(say::DEBUG);

    fastNLOAlphas fnloOld("theorFiles/fastnlo-cms-incjets-arxiv-1605.04436-xsec001.tab", "CT14nlo", 0);
    fastNLOAlphas fnloNew("theorFiles/InclusiveNJets_fnl5362h_v23_fix.tab", "CT14nlo", 0);
    fnloOld.SetContributionON(fastNLO::kFixedOrder,0,true);
    fnloOld.SetContributionON(fastNLO::kFixedOrder,1,true);
    fnloOld.SetUnits(fastNLO::kPublicationUnits);

    fnloNew.SetContributionON(fastNLO::kFixedOrder,0,true);
    fnloNew.SetContributionON(fastNLO::kFixedOrder,1,true);
    fnloNew.SetUnits(fastNLO::kPublicationUnits);



    //vector<vector<TH1D*>> histNewScl = getScaleuncHistos(fnloNew);
    //vector<vector<TH1D*>> histNewPDF = getPDFuncHistos(fnloNew);
    //vector<vector<TH1D*>> histOldScl = getScaleuncHistos(fnloOld);
    //vector<vector<TH1D*>> histOldPDF = getPDFuncHistos(fnloOld);


    //SaveHistos(histNewPDF, "histNewCT14PDF");
    //SaveHistos(histNewScl, "histNewCT14Scl");
    //SaveHistos(histOldPDF, "histOldCT14PDF");
    //SaveHistos(histOldScl, "histOldCT14Scl");

    //histNewScl[0][0]->Print("all");
    //histNewScl[0][0]->Write();
    //fOut->Write();
    //fOut->Close();

    /*
    return 0;

    int nMem = fnloOld.GetNPDFMembers();


    vector<TH1D*> oldH = readHisto(fnloOld);
    vector<TH1D*> newH = readHisto(fnloNew);


    TCanvas *can = new TCanvas(rn(), "");
    gStyle->SetOptStat(0);

    DividePad( {1,1,1,1,1}, {1});

    for(int yB = 0; yB < 5; ++yB) {

        TH1D *oldHH = rebin(oldH[yB], newH[yB]);


        oldHH->Divide(newH[yB]);

        can->cd(yB+1);
        gPad->SetLogx();
        oldHH->Draw();
        GetYaxis()->SetRangeUser(0.99, 1.01);
        GetXaxis()->SetTitle("p_{T} [GeV]");
        GetXaxis()->SetRangeUser(100, 2600);

        for(int i = 1; i < oldHH->GetNbinsX(); ++i) {
            cout << i << " "<<  oldHH->GetBinLowEdge(i) <<" "<<  oldHH->GetBinContent(i) <<" "<< newH[yB]->GetBinContent(i) <<  endl;
        }

    }
    can->SaveAs("tableCheck.pdf");
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
