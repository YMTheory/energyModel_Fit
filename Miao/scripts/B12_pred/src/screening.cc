#include "global.hh"
#include "screening.hh"

#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>

screening* gScreen = (screening*)0;

screening::screening()
{
    mScreen = new TF1("Screen", "", 0, g.E0);
}

double screening::screenCorr(double betaE) {

    double betaW = betaE/0.511+1;
    double betaP = TMath::Sqrt(betaW*betaW-1);
    double Z_bar = g.Z;
    double V0 = g.alpha*g.alpha*TMath::Power(Z_bar, 0.75)*1.4584000;  // linear interpolation from table
    double betaW_bar = betaW - V0;
    double betaP_bar = TMath::Sqrt(betaW_bar*betaW_bar-1);
    double y = g.alpha*g.Z*betaW/betaP;
    double y_bar = g.alpha*g.Z*betaW_bar/betaP_bar;
    double S1 = betaW_bar/betaW * TMath::Power(betaP_bar/betaP, 2*g.gamma-1)*TMath::Exp(TMath::Pi()*(y_bar-y)) * 0.253564;
    double Z0K0 = 1;
    double Z2K0 = 0.577216; double Z2K2 = -1.644934;
    double Z4K0 = 0.722126; double Z4K2 = -2.151539; double Z4K4 = 1.894066;
    double Z6K0 = 0.730658; double Z6K2 = -2.993953; double Z6K4 = 4.107516; double Z6K6 = -1.971102;
    double tt1 = g.alpha*g.Z; double tt2 = betaW_bar/betaP_bar;
    double S2 = Z0K0*TMath::Power(tt1, 0)*TMath::Power(tt2, 0)
        + Z2K0*TMath::Power(tt1,2)*TMath::Power(tt2,0)
        + Z2K2*TMath::Power(tt1,2)*TMath::Power(tt2,2)
        + Z4K0*TMath::Power(tt1,4)*TMath::Power(tt2,0)
        + Z4K2*TMath::Power(tt1,4)*TMath::Power(tt2,2)
        + Z4K4*TMath::Power(tt1,4)*TMath::Power(tt2,4)
        + Z6K0*TMath::Power(tt1,6)*TMath::Power(tt2,0)
        + Z6K2*TMath::Power(tt1,6)*TMath::Power(tt2,2)
        + Z6K4*TMath::Power(tt1,6)*TMath::Power(tt2,4)
        + Z6K6*TMath::Power(tt1,6)*TMath::Power(tt2,6);
    double S = 1;
    if(betaW >= V0) { S = S1+S2; }

    return S;
}

void screening::Plot() {

    TGraph* graphScreen = new TGraph();
    const int nPoints = 1000;
    for(int iPoint = 0; iPoint<nPoints; iPoint++) {
        double E = g.E0/nPoints*(iPoint+1);
        graphScreen->SetPoint(iPoint, E, screenCorr(E));
    }
    graphScreen->SetLineColor(kBlue+1);
    graphScreen->SetName("Screening");

    TFile* file = new TFile("screening.root", "RECREATE");
    graphScreen->Write();
    file->Close();
}

double gScreenFcn(double* x, double* p) {
    double betaE = x[0];
    if(gScreen == 0) gScreen = new screening();
    return gScreen->screenCorr(betaE);
}
