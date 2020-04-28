#include "junoParameters.hh"

typedef junoParameters JUNOP;
typedef ScintillatorParameterization SP;

SP JUNOP::scintillatorParameterization = kSimulation;

std::string JUNOP::stopPow_File       = "data/electron/StopPow1.txt";
std::string JUNOP::quenchNL_File      = "data/electron/Quench.root";
std::string JUNOP::cerenkovNL_File    = "data/electron/CerenkovPE.txt";
std::string JUNOP::electronLSNL_File  = "data/electron/electron_total.txt";
std::string JUNOP::electronOut_File   = "output/electron/electronFit.root";

std::string JUNOP::gammaLSNL_File     = "data/naked_gamma/gamma.dat";
std::string JUNOP::gammaPdf_File      = "data/Gamma_Electron.root";
std::string JUNOP::gammaOut_File      = "output/gamma/gammaFit.root";
double JUNOP::m_gammaError            = 1.;
