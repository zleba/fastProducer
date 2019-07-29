///********************************************************************
///
///     fastNLO_toolkit: fnlo-tk-example
///     Example program for the user to test cross-section calculations
///     with fastNLO.
///
///     D. Britzger, K. Rabbertz
///
///********************************************************************

#include "config.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "fastnlotk/fastNLOLHAPDF.h"

//! Includes for filling ROOT histograms
//! Usable only when configured with '--with-root=/path/to/root' option
#ifdef WITH_ROOT
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#endif
//! End of ROOT part

//! Function prototype for flexible-scale function
double Function_Mu(double s1, double s2);

//______________________________________________________________________________________________________________
int main(int argc, char** argv) {

   //! --- namespaces
   using namespace std;
   using namespace say;       //! namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;   //! namespace for fastNLO constants

   //! --- Set verbosity level
   SetGlobalVerbosity(INFO);

   //! --- Print program purpose
   yell << _CSEPSC << endl;
   info["fnlo-tk-example"] << "Example program for the user to test cross-section calculations" << endl;
   info["fnlo-tk-example"] << "with fastNLO" << endl;
   yell << _SSEPSC << endl;
   info["fnlo-tk-example"] << "For more explanations type:" << endl;
   info["fnlo-tk-example"] << "./fnlo-tk-example -h" << endl;
   info["fnlo-tk-example"] << "and consult the provided source code 'fnlo-tk-example.cc'." << endl;
   yell << _CSEPSC << endl;

   //! ---  Parse commmand line
   yell << "" << endl;
   yell << _CSEPSC << endl;
   shout["fnlo-tk-example"] << "fastNLO Example Evaluator" << endl;
   yell << _SSEPSC << endl;
   string tablename;
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   string PDFFile = "CT10nlo";
#else
   string PDFFile = "CT10nlo.LHgrid";
#endif
   if (argc <= 1) {
      error["fnlo-tk-example"] << "No fastNLO table specified!" << endl;
      shout["fnlo-tk-example"] << "For an explanation of command line arguments type:" << endl;
      shout["fnlo-tk-example"] << "./fnlo-tk-example -h" << endl;
      yell << _CSEPSC << endl;
      exit(1);
   } else {
      //! --- fastNLO table
      tablename = (const char*) argv[1];
      //! --- Usage info
      if (tablename == "-h") {
         yell << " #" << endl;
         info["fnlo-tk-example"] << "This is an example program to evaluate a fastNLO table" << endl;
         info["fnlo-tk-example"] << "that a user can check and adapt to his liking" << endl;
         man << "" << endl;
         man << "Usage: ./fnlo-tk-example <fastNLOtable.tab> [PDF]" << endl;
         man << "       Arguments: <> mandatory; [] optional." << endl;
         man << "<fastNLOtable.tab>: Table input file fnl2342b.tab" << endl;
         man << "[PDF]: PDF set, def. = CT10nlo" << endl;
         man << "   For LHAPDF5: Specify set names WITH filename extension, e.g. \".LHgrid\"." << endl;
         man << "   For LHAPDF6: Specify set names WITHOUT filename extension." << endl;
         man << "   If the PDF set still is not found, then:" << endl;
         man << "   - Check, whether the LHAPDF environment variable is set correctly." << endl;
         man << "   - Specify the PDF set including the absolute path." << endl;
         man << "   - Download the desired PDF set from the LHAPDF web site." << endl;
         yell << " #" << endl;
         yell << _CSEPSC << endl;
         return 0;
      } else {
         shout["fnlo-tk-example"] << "Evaluating table: "  <<  tablename << endl;
      }
   }
   //! --- PDF choice
   if (argc > 2) {
      PDFFile = (const char*) argv[2];
   }

   //! --- Give some output
   shout["fnlo-tk-example"] << "Evaluating table: " << tablename << endl;
   shout["fnlo-tk-example"] << "Using PDF set   : " << PDFFile << endl;

   //! --- This is your playgroud to use fastNLO
   //!     Calculate cross setions and/or test some options
   //!     For some explanation and function calls, please see
   //!     the other code examples in './src/' and
   //!     the Doxygen documentation.

   //! --- Example calculation
   fastNLOLHAPDF fnlo(tablename,PDFFile,0);     //! initialize a fastNLO instance with interface to LHAPDF.
   fnlo.PrintContributionSummary(0);             //! print some valuable information
   //fnlo.Print();                              //! print even more information
   //fnlo.SetUnits(kAbsoluteUnits);             //! Use units as specified in the publication or in barns.
   //fnlo.SetContributionON(fastNLO::kFixedOrder,0,false); //! switch contributions on/off. By default LO and NLO.
   //fnlo.SetContributionON(fastNLO::kFixedOrder,1,true);
   //fnlo.SetContributionON(fastNLO::kFixedOrder,2,true); //! NNLO must be switched on explicitly
   fnlo.CalcCrossSection();                     //! Calculate the cross section
   fnlo.PrintCrossSections();                   //! Print cross section to screen

   vector<double> xs = fnlo.GetCrossSection();  //! Access cross sections for later usage

   //! Finish?
   //return 0;


   //! --- Example calculation of cross section including relative uncertainty
   //!
   //! Enumerators for type of scale and PDF uncertainty
   //! Choices for scale uncertainty are:
   //!   kScaleNone           : no scale uncertainty, only central scale (mu_r,mu_f) = (1,1) evaluated
   //!   kSymmetricTwoPoint   : symmetric (mu_r,mu_f) scale variations by factors (1/2,1/2), (2,2)
   //!   kAsymmetricSixPoint  : asymmetric (mu_r,mu_f) scale variations by factors (1/2,1/2), (2,2) plus
   //!                          (1/2,1), (1,1/2), (1,2), (2,1)
   //! Choices for PDF uncertainty are:
   //!   kPDFNone             : No PDF uncertainty, only averaged cross section result evaluated (Correct for NNPDF, wrong otherwise!)
   //!   kLHAPDF6             : LHAPDF6 uncertainties (recommended if LHAPDF6 is available)
   //!                          Uses standard formula for NNPDF with symmetric uncertainties around mean.
   //!                          For alternative see specific calls below and the LHAPDF6 documentation.
   //!   kHessianSymmetric    : symmetric Hessian PDF uncertainties (ABM, (G)JR)
   //!   kHessianAsymmetric   : asymmetric Hessian PDF uncertainties
   //!   kHessianAsymmetricMax: asymmetric Hessian PDF uncertainties with pairwise max deviations per eigenvector (CTEQ,MRST|MSTW)
   //!   kHessianCTEQCL68     : like kHessianAsymmetricMax, but with uncertainties rescaled to CL68
   //!   kMCSampling          : statistical sampling PDF uncertainties and central value defined as mean (NNPDF)
   //!   kHeraPDF10           : HERAPDF 1.0 uncertainties (NOT implemented yet, neither here nor in LHAPDF6!)

   EScaleUncertaintyStyle eScaleUnc = kAsymmetricSixPoint;
   //   EPDFUncertaintyStyle   ePDFUnc   = kHessianCTEQCL68;
   // Return values are three vectors xs, dxsu, dxsl in struct XsUnc
   XsUncertainty XsUnc;
   XsUnc = fnlo.GetScaleUncertainty(eScaleUnc);
   //   XsUnc  = fnlo.GetPDFUncertainty(ePDFUnc);

   cout << _CSEPSC << endl;
   cout << " # Relative Scale Uncertainties (6P)" << endl;
   cout << " # bin      cross section           lower uncertainty       upper uncertainty" << endl;
   cout << _SSEPSC << endl;
   for ( unsigned int iobs=0;iobs<XsUnc.xs.size();iobs++ ) {
      printf("%5.i      %#18.11E      %#18.11E      %#18.11E\n",iobs+1,XsUnc.xs[iobs],XsUnc.dxsl[iobs],XsUnc.dxsu[iobs]);
   }
   cout << _CSEPSC << endl;

   //! Finish?
   return 0;


   //! --- Example calculation of cross section including PDF uncertainty using LHAPDF6 specifically
   //!
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   vector<LHAPDF::PDFUncertainty> PDFUnc = fnlo.GetPDFUncertaintyLHAPDF();
   //! Using the following call the CL can be changed and the alternativ asymmetric NNPDF uncertainty
   //! around the median as implemented in LHAPDF6 can be switched on.
   //   vector<LHAPDF::PDFUncertainty> PDFUnc = fnlo.GetPDFUncertaintyLHAPDF(100*erf(1/sqrt(2)),true);
   vector<double> errup    = fnlo.CalcPDFUncertaintyRelPlus(PDFUnc);
   vector<double> errdn    = fnlo.CalcPDFUncertaintyRelMinus(PDFUnc);
   vector<double> errupabs = fnlo.CalcPDFUncertaintyPlus(PDFUnc);
   vector<double> errdnabs = fnlo.CalcPDFUncertaintyMinus(PDFUnc);
   vector<double> central  = fnlo.CalcPDFUncertaintyCentral(PDFUnc);

   cout << _CSEPSC << endl;
   cout << " # Relative and Absolute PDF Uncertainties via LHAPDF6" << endl;
   cout << " # bin      cross section      lower rel. uncertainty   upper rel. uncertainty   lower abs. uncertainty   upper abs. uncertainty" << endl;
   cout << _SSEPSC << endl;
   for ( unsigned int iobs=0;iobs<central.size();iobs++ ) {
      printf("%5.i      %#18.11E      %#18.11E      %#18.11E      %#18.11E      %#18.11E\n",iobs+1,central[iobs],errdn[iobs],errup[iobs],errdnabs[iobs],errupabs[iobs]);
   }
   cout << _CSEPSC << endl;

   //! Finish?
   return 0;
#endif


   //! --- Example filling of ROOT histogram with previously calculated cross section and uncertainty
   //! Usable only when configured with '--with-root=/path/to/root' option
#ifdef WITH_ROOT
   TString out_file_name = "./fnlo_out.root";
   TFile *file_out = new TFile(out_file_name,"NEW");
   TH1D *histo1 = new TH1D("Cross Section Bins","fastNLO",(int)xs.size()+1,0.5,xs.size()+0.5);
   histo1->GetXaxis()->SetTitle("Bin Number");
   histo1->GetYaxis()->SetTitle("Cross Section");
   for( unsigned int iobs=0;iobs<xs.size();iobs++ ){
      histo1->SetBinContent(iobs+1,xs[iobs]);
      // Symmetrize uncertainty since ROOT does not support histograms with asymmetric errors
      histo1->SetBinError(iobs+1,sqrt(XsUnc.dxsl[iobs]*XsUnc.dxsl[iobs] + XsUnc.dxsu[iobs]*XsUnc.dxsu[iobs])*xs[iobs]/2);
   }

   file_out->cd();
   file_out->Write();
   file_out->Close();
   //! End of ROOT part

   //! Finish?
   return 0;
#endif


   //! Example code how to loop over all PDF eigenvectors
   cout<<"\n fnlo-tk-example: Now we want to loop over the eigenvectors of "<<PDFFile<<"."<<endl<<endl;
   fnlo.SetLHAPDFFilename(PDFFile); //! we use again the 'nominal' PDF-file
   int nEig = fnlo.GetNPDFMembers(); //! How many eigenvectors are there?
   cout<<" fnlo-tk-example: There are "<<nEig<<" Eigenvalue sets in "<<PDFFile<<endl;
   for ( int i = 0 ; i<nEig ; i++ ) { //! start with 0
      cout<<" fnlo-tk-example: Setting PDF member: "<<i<<" ***"<<endl;
      fnlo.SetLHAPDFMember(i);  //! specify the PDF member
      fnlo.CalcCrossSection();  //! redo cross section calculation
      fnlo.PrintCrossSections(); //! print new cross sections to screen
      //! write tot file
      vector<double> cs = fnlo.GetCrossSection(); //! get cross setions for further usage (i.e. printing to file)
   }

   //! Finish
   return 0;

   //! Example to use different scale factors
   cout<<"\n fnlo-tk-example: Now we use a scale factor of 2."<<endl;
   fnlo.SetScaleFactorsMuRMuF(2.0,2.0);
   fnlo.CalcCrossSection();
   fnlo.PrintCrossSections();


   cout<<"\n fnlo-tk-example: Now we use a scale factor of 0.5."<<endl;
   fnlo.SetScaleFactorsMuRMuF(0.5,0.5);
   fnlo.CalcCrossSection();
   fnlo.PrintCrossSections();


   cout<<"\n fnlo-tk-example: Now we go back to the nominal result: 1."<<endl;
   fnlo.SetScaleFactorsMuRMuF(1.0,1.0);
   fnlo.CalcCrossSection();
   fnlo.PrintCrossSections();

   //! Finish
   return 0;

}


double Function_Mu(double s1, double s2) {
   //! --- fastNLO user: This is an example function
   //!     to demonstrate how you might perform the
   //!     definition of the scales using a
   //!     'flexible-scale'-table, where a function
   //!     of s1 and s2 can be used.
   //!     Which variables s1 and s2 stand for are
   //!     coded in the fastNLO table.
   double mu = 173.;
   return mu;
}
