// Author: Daniel Britzger
// DESY, 20/04/2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO_reader_2.1.0                                                //
//  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch           //
//                                                                      //
//  The projects web page can be found at:                              //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//  If you use this code, please cite:                                  //
//    T. Kluge, K. Rabbertz and M. Wobisch, hep-ph/0609285              //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch,        //
//       arXiv:1109.1310                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//
//  fastNLOAlphas
//  This class inherits the PDF interface from
//  fastNLOLHAPDF, while the alpha_s evolution
//  is superseeded by the Alphas.h class.
//lhasub
//////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
//#include "fastnlotk/fastNLOReader.h"
//#include "fastnlotk/speaker.h"
#include "fastnlotk/fastNLOLHAPDF.h"
#include "fastnlotk/fastNLOHoppet.h"
#include "fastnlotk/HoppetInterface.h"

using namespace std;

//______________________________________________________________________________
//
//
fastNLOHoppet::fastNLOHoppet(std::string name) : fastNLOLHAPDF(name) {
   //Set some meaningful initial values
   SetPDGValues();
   // KR: Note: LHAPDF values cannot be taken here, since the class instantiation may
   //     happen before defining the PDF set!
   //   SetLHAPDFValues();
};

fastNLOHoppet::fastNLOHoppet(std::string name, std::string LHAPDFFile, int PDFSet = 0) : fastNLOLHAPDF(name,LHAPDFFile,PDFSet) {
   //Set some meaningful initial values
   SetPDGValues();
   // KR: For consistency with usage above.
   //   SetLHAPDFValues();
};



// Getters
double fastNLOHoppet::GetQMass(int pdgid) const {
   return HoppetInterface::QMass[pdgid];
}
double fastNLOHoppet::GetMz() const {
   return HoppetInterface::fMz;
}
int fastNLOHoppet::GetNFlavor() const {
   return HoppetInterface::fnFlavor;
}
int fastNLOHoppet::GetNLoop() const {
   return HoppetInterface::fnLoop;
}
double fastNLOHoppet::GetAlphasMz() const {
   return HoppetInterface::fAlphasMz;
};



// Setters
void fastNLOHoppet::SetQMass(int pdgid, double qmass) {
   HoppetInterface::QMass[pdgid] = qmass;
   HoppetInterface::InitHoppet(*this);
}
void fastNLOHoppet::SetMz(double Mz) {
   HoppetInterface::fMz = Mz;
   HoppetInterface::InitHoppet(*this);
}
void fastNLOHoppet::SetNFlavor(int nflavor) {
   HoppetInterface::fnFlavor = nflavor;
   HoppetInterface::InitHoppet(*this);
}
void fastNLOHoppet::SetNLoop(int  nloop) {
   HoppetInterface::fnLoop = nloop;
   HoppetInterface::InitHoppet(*this);
}
void fastNLOHoppet::SetAlphasMz(double AlphasMz, bool ReCalcCrossSection) {
   HoppetInterface::fAlphasMz    = AlphasMz;
   HoppetInterface::InitHoppet(*this);
}



// Combined Setters
void fastNLOHoppet::SetPDGValues() {
   // Initialize with PDG values
   HoppetInterface::QMass[0]  = PDG_MD;
   HoppetInterface::QMass[1]  = PDG_MU;
   HoppetInterface::QMass[2]  = PDG_MS;
   HoppetInterface::QMass[3]  = PDG_MC;
   HoppetInterface::QMass[4]  = PDG_MB;
   HoppetInterface::QMass[5]  = PDG_MT;
   HoppetInterface::fMz       = PDG_MZ;
   //Variable flavor number scheme
   HoppetInterface::fnFlavor = 0;
   //2-loop alpha_s evolution
   HoppetInterface::fnLoop = 2;
   HoppetInterface::fAlphasMz = PDG_ASMZ;
   HoppetInterface::InitHoppet(*this);
}

void fastNLOHoppet::SetLHAPDFValues() {
   //Be sure LHAPDF is initialized when reading the properties
   if (fchksum == 0 || fchksum != CalcChecksum(1.)) {
      if ( ! InitPDF() ) {
         logger.error["SetLHAPDFValues"]<<"No LHAPDF set initialized, aborting!\n";
         exit(1);
      } else {
         FillPDFCache();
      }
   }

   for (int i = 0; i < 6; i++) {
      HoppetInterface::QMass[i] = LHAPDF::getQMass(i+1);
   }
   //How to read LHAPDF Mz???
   HoppetInterface::fMz = PDG_MZ;
   HoppetInterface::fnFlavor = LHAPDF::getNf();
   HoppetInterface::fnLoop = LHAPDF::getOrderAlphaS();
   HoppetInterface::fAlphasMz = LHAPDF::alphasPDF(HoppetInterface::fMz);
   HoppetInterface::InitHoppet(*this);
}

bool fastNLOHoppet::InitPDF() {
   bool init = fastNLOLHAPDF::InitPDF();
   HoppetInterface::InitHoppet(*this);
   return init;
}

// Evolution
double fastNLOHoppet::EvolveAlphas(double Q ) const {
   return HoppetInterface::EvolveAlphas(Q);
}

std::vector<double> fastNLOHoppet::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return HoppetInterface::GetXFX(xp, muf);
}
