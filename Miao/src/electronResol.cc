#include "electronResol.hh"

#include <iostream>

#include "TMath.h"
#include "TRandom.h"

using namespace std;

double electronResol::m_pA = 2.57576e-02;
double electronResol::m_pB = 6.86431e-03;
double electronResol::m_pC = 0.;

double electronResol::energySmear(double Evis) {
    double resol = TMath::Sqrt(m_pA*m_pA/Evis+m_pB*m_pB+m_pC*m_pC/Evis/Evis);
    double Esmear = gRandom->Gaus(Evis, Evis*resol);
    return Esmear;
}

double electronResol::Resolution(double Evis) {
    double resol = TMath::Sqrt(m_pA*m_pA/Evis+m_pB*m_pB+m_pC*m_pC/Evis/Evis);
    return resol;
}
