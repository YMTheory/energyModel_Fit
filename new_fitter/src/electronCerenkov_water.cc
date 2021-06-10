#include "electronCerenkov_water.hh"

#include <TMath.h>

double electronCerenkov_water::m_w1 = 0;
double electronCerenkov_water::m_w2 = 0;
double electronCerenkov_water::m_w3 = 0;
double electronCerenkov_water::m_w4 = 0;

double electronCerenkov_water::m_r0 = 0;
double electronCerenkov_water::m_r1 = 0;
double electronCerenkov_water::m_r2 = 0;

double gCerPE_water(double* x, double* p)
{
    double E = x[0];
    double w1 = p[0];
    double w2 = p[1];
    double w3 = p[2];
    double w4 = p[3];
    double xx = TMath::Log(1+E/0.2);
    double npe = (w1*xx + w2*xx*xx + w3*xx*xx*xx) * (1/E + w4) * E   ;
    return npe;

}

double gCerPESigma_water(double* x, double* p)
{
    double N = x[0];
    double r0 = p[0];
    double r1 = p[1];
    double r2 = p[2];

    double sigma2 = r0 + r1*N + r2*N*N;
    return TMath::Sqrt(sigma2);
}

TF1* electronCerenkov_water::fCerPE_water = new TF1("fCerPE_water", gCerPE_water, 0, 60, 4);
TF1* electronCerenkov_water::fCerPESigma_water = new TF1("fCerPESigma_water", gCerPESigma_water, 0, 1200, 3);



electronCerenkov_water::electronCerenkov_water()
{;}

electronCerenkov_water::~electronCerenkov_water()
{;}


double electronCerenkov_water::getCerenkovPE(double E)
{
    return fCerPE_water->Eval(E);
}



double electronCerenkov_water::getCerenkovPESigma(double N) 
{
    return fCerPESigma_water->Eval(N);
}
