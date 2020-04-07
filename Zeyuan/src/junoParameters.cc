#include "junoParameters.hh"

typedef junoParameters JUNOP;

int JUNOP::nFitParameter = 3;
int JUNOP::fitPrintLevel = 1;

// MC nominal
double JUNOP::p0_start       = 0.98;
double JUNOP::p1_start       = 6.5e-3;
double JUNOP::p2_start       = 1 ;

string JUNOP::gammaPdf_file   = "data/Gamma_Electron.root";
string JUNOP::gammaData_file   = "data/naked_gamma/gamma.dat";

bool JUNOP::fitGamma  = true;

double JUNOP::gamScale       = 1.0;
double JUNOP::gamScaleError  = 0.0015;

bool JUNOP::fixScintP0   = false; //energy scale
bool JUNOP::fixScintP1   = false; //kb
bool JUNOP::fixScintP2   = false; //Cherenkov
bool JUNOP::fixGamScale  = true;

