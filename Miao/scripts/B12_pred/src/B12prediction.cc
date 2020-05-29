#include "global.hh"
#include "B12prediction.hh"
#include "fermiFunc.hh"
#include "finite_size.hh"
#include "screening.hh"
#include "weak_magnetism.hh"

#include <TGraph.h>
#include <TFile.h>

using namespace std;

double B12prediction::predSpec(double betaE) {
    fermiFunc* fermi = new fermiFunc();
    finite_size* fs = new finite_size();
    screening* screen = new screening();
    weak_magnetism* wm = new weak_magnetism();

    return fermi->phasespace(betaE)*fermi->fermiCorr(betaE)*fs->elecmagCorr(betaE)*fs->weakCorr(betaE)*screen->screenCorr(betaE)*wm->weakMagCorr(betaE);
}

void B12prediction::Plot() {
    TGraph* graphSpec = new TGraph();
    const int nPoints = 200;
    double maxPoint = -1;
    for(int iPoint=0; iPoint<nPoints; iPoint++) {
        double E = g.E0/nPoints*(iPoint+1);
        if(predSpec(E)>maxPoint) maxPoint = predSpec(E);
    }
    
    cout << "Maximum: " << maxPoint << endl;
    for(int iPoint=0; iPoint<nPoints; iPoint++) {
        double E = g.E0/nPoints*(iPoint+1);
        graphSpec->SetPoint(iPoint, E, predSpec(E)/maxPoint);
    }

    graphSpec->SetMarkerColor(kBlue);
    graphSpec->SetMarkerStyle(20);
    graphSpec->SetMarkerSize(0.2);
    graphSpec->SetName("predSpec");
    
    TFile* file = new TFile("predSpec.root", "recreate");
    graphSpec->Write();
    file->Close();
}
