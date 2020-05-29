#include "junoParameters.hh"

typedef junoParameters JUNOP;
typedef ScintillatorParameterization SP;

SP JUNOP::scintillatorParameterization = kSimulation; /*kEmpirical;*/ /*kIntegral;*/

std::string JUNOP::stopPow_File       = "data/electron/StopPow1.txt";
std::string JUNOP::quenchNL_File      = "data/electron/Quench1.root";
std::string JUNOP::cerenkovNL_File    = "data/electron/Cer.dat";
std::string JUNOP::quenchNL_outFile   = "output/electron/QuenchNL.root";
std::string JUNOP::cerenkov_outFile   = "output/electron/CerenkovNL.root";
//std::string JUNOP::cerenkovNL_File    = "data/electron/CerenkovPE.txt";
std::string JUNOP::electronLSNL_File  = "data/electron/electron_total.txt";
std::string JUNOP::electronOut_File   = "output/electron/electronFit.root";

std::string JUNOP::gammaLSNL_File     = "data/naked_gamma/gamma.dat";
std::string JUNOP::gammaPdf_File      = "data/Gamma_Electron.root";
std::string JUNOP::gammaOut_File      = "output/gamma/gammaFit.root";
double JUNOP::m_gammaError            = 1.;

std::string JUNOP::B12DataFile        = "data/electron/B12.root";
std::string JUNOP::B12PredFile        = "data/electron/B12_beta.txt";   
std::string JUNOP::B12CalcFile        = "output/electron/B12_Calc.root";

std::string JUNOP::B12_File           = "data/electron/B12_file.root";
std::string JUNOP::B12Out_File        = "output/electron/B12Fit.root";
std::string JUNOP::B12Spec_File       = "data/electron/B12_nonl.txt";

std::string JUNOP::C11_File           = "data/electron/C11_file.root";
std::string JUNOP::C11Out_File        = "output/electron/C11Fit.root";
std::string JUNOP::C11Spec_File       = "data/electron/C11_nonl.txt";
