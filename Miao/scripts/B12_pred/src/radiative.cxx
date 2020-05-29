#include "global.hh"
#include "radiative.hh"

#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>

radiative* gRadiative = (radiative*)0;

radiative::radiative() {
    mRadiative = new TF1("radiative", "", 0, g.E0);
}

double radiative::radiaCorr(double betaE) {
    
}

void radiative::Plot()
{}


double gRadiaFcn(double* x, double* p) {
    double betaE = x[0];
    if (gRadiative == 0) gRadiative = new radiative();
    return gRadiative->radiaCorr(betaE);
}
