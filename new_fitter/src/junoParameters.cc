#include "junoParameters.hh"

typedef junoParameters JUNOP;
typedef ScintillatorParameterization SP;

SP JUNOP::scintillatorParameterization = kSimulation;  /*kEmpirical;*/ /*kIntegral;*/

std::string JUNOP::stopPow_File       = "data/electron/StopPow1.txt";
std::string JUNOP::quenchNL_File      = "data/electron/Quench5.root";
std::string JUNOP::cerenkovNL_File    = "data/electron/Cer.dat";
std::string JUNOP::quenchNL_outFile   = "output/electron/QuenchNL.root";
std::string JUNOP::cerenkov_outFile   = "output/electron/CerenkovNL.root";
//std::string JUNOP::cerenkovNL_File    = "data/electron/CerenkovPE.txt";
std::string JUNOP::electronLSNL_File  = "data/electron/electron_total.txt";
std::string JUNOP::electronOut_File   = "output/electron/electronFit.root";

std::string JUNOP::m_calcOption       = "twolayer";
std::string JUNOP::m_nonlMode         = "analytic";

std::string JUNOP::gammaLSNL_File     = "data/gamma/gamma1.dat";
std::string JUNOP::gammaPdf_File      = "data/gamma/Gamma_Electron1.root";
std::string JUNOP::gammaOut_File      = "output/gamma/gammaFit.root";
double JUNOP::m_gammaError            = 1.;
bool JUNOP::fitGammaSources           = true;
std::string JUNOP::gammaNLOption      = "prmelec";  

bool JUNOP::fitB12                    = false;
std::string JUNOP::B12DataFile        = "data/electron/B12.root";
std::string JUNOP::B12PredFile        = "data/electron/B12_beta.txt";   
std::string JUNOP::B12CalcFile        = "output/electron/B12_Calc.root";
double JUNOP::b12FitMinE              = 4;
double JUNOP::b12FitMaxE              = 14;
double JUNOP::b12VertexCut            = 0.80;

bool JUNOP::fitC11                    = false;
std::string JUNOP::B12_File           = "data/electron/B12_file.root";
std::string JUNOP::B12Out_File        = "output/electron/B12Fit.root";
std::string JUNOP::B12Spec_File       = "data/electron/B12_nonl.txt";

std::string JUNOP::C11DataFile        = "data/electron/C11.root";
std::string JUNOP::C11_File           = "data/electron/C11_file.root";
std::string JUNOP::C11Out_File        = "output/electron/C11Fit.root";
std::string JUNOP::C11Spec_File       = "data/electron/C11_nonl.txt";
double JUNOP::c11FitMinE              = 0.9;
double JUNOP::c11FitMaxE              = 1.8;

bool JUNOP::fitC10                    = false;
std::string JUNOP::C10DataFile        = "data/electron/C10.root";
std::string JUNOP::C10Out_File        = "output/electron/C10Fit.root";
double JUNOP::c10FitMinE              = 1.6;
double JUNOP::c10FitMaxE              = 3.4;