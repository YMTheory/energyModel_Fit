#include "electronQuench.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <TMath.h>
#include <TGraph.h>

using namespace std;

double electronQuench::m_birk1     = 1.;

double electronQuench::s_quenchingShape1[s_nKb][s_nSamples] = {0};

electronQuench::electronQuench()
{;}

electronQuench::~electronQuench()
{;}

double electronQuench::Integral_BirkLaw(double kB1, double E)
{
    if( Etrue.size()==0 || StopPow.size() == 0) {
        std::cout << " >>> No Data Loaded in Vector <<< " << endl;
    } else if (Etrue.size() != StopPow.size() ) {
        std::cout << " >>> Quench Integral Vectors have Different Lengths <<< " << std::endl;
    } else {
        
    }
    
    return 1.0;
}




