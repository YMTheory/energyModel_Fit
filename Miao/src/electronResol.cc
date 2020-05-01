#include "electronResol.hh"

#include <iostream>

#include "TMath.h"
#include "TRandom.h"

using namespace std;

double electronResol::m_pA = 0.026;
double electronResol::m_pB = 0.0068;
double electronResol::m_pC = 0.;

double electronResol::energySmear(double Evis) {
    double resol = TMath::Sqrt(m_pA*m_pA/Evis+m_pB*m_pB+m_pC*m_pC/Evis/Evis);
    double Esmear = gRandom->Gaus(Evis, Evis*resol);
    return Esmear;
}
