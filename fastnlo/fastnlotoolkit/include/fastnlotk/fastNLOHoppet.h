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

#ifndef FASTNLOHOPPET
#define FASTNLOHOPPET

//#include "fastNLOReader.h"
//#include <iostream>
//#include <cstdio>
//#include <cstdlib>
//#include <LHAPDF/LHAPDF.h>
//#include "speaker.h"
#include "fastNLOLHAPDF.h"
//#include "hoppet_v1.h"



class fastNLOHoppet : public fastNLOLHAPDF {

   public:
      fastNLOHoppet(std::string name);
      fastNLOHoppet(std::string name, std::string LHAPDFFile, int PDFSet);
      // ---- Alphas vars ---- //
      // Setters
      void SetMz(double Mz);
      void SetNFlavor(int nflavor);
      void SetNLoop(int nloop);
      void SetQMass(int pdgid, double qmass);
      void SetAlphasMz(double AlphasMz, bool ReCalcCrossSection = false);
      void SetLHAPDFValues();
      void SetPDGValues();
      virtual bool InitPDF();
      // Getters
      double GetMz() const;
      double GetQMass(int pdgid) const;
      int GetNFlavor() const;
      int GetNLoop() const;
      double GetAlphasMz() const;



   protected:

      // inherited functions
      virtual double EvolveAlphas(double Q) const ;
      //bool InitPDF();
      virtual std::vector<double> GetXFX(double xp, double muf) const ;
      // ---- Alphas vars ---- //
};

#endif
