#include "global.hh"
#include "weak_magnetism.hh"

#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>

weak_magnetism::weak_magnetism()
{
    mWeakMag = new TF1("WeakMagnetism", "", 0, g.E0);
}

double weak_magnetism::weakMagCorr(double betaE) {
    double betaW = betaE/0.511+1;
    double delta = 0.0048;
    return 1+delta*betaW;

}

void weak_magnetism::Plot()
{
    TGraph* graphWM = new TGraph();
    const int nPoints = 1000;
    for(int iPoint=0; iPoint<nPoints; iPoint++) {
        double E = g.E0/nPoints*(iPoint+1);
        graphWM->SetPoint(iPoint, E, weakMagCorr(E));
    }
    graphWM->SetName("WeakMagnetism");
    graphWM->SetLineColor(kRed+1);

    TFile *file = new TFile("WeakMagnetism.root", "RECREATE");
    graphWM->Write();
    file->Close();
}
