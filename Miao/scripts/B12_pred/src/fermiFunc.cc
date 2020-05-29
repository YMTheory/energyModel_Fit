#include "global.hh"
#include "fermiFunc.hh"

#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>

fermiFunc* gFermiFunc  = (fermiFunc*)0;

fermiFunc::fermiFunc()
{
    mFermiFcn = new TF1("fermi_fun", gFermiFcn, 0, g.E0);
}

double fermiFunc::phasespace(double betaE) {

    // traditional phase space : 
    double betaW = betaE/0.511+1;
    double betaP = TMath::Sqrt(betaW*betaW-1);
    double sp = betaP*betaP*(g.W0-betaW)*(g.W0-betaW);
    return sp;
}

double fermiFunc::fermiCorr(double betaE) {

    double betaW = betaE/0.511+1;
    double betaP = TMath::Sqrt(betaW*betaW-1);

    double F1 = 2*(g.gamma+1)*TMath::Power(2*g.R, 2*(g.gamma-1))*0.253564; // energy-independent part
    // Table for expansions ...
    double Z0K0L0 = 1; 
    double Z1K1L0 = 3.141593;
    double Z2K0L0 = 0.577216; double Z2K0L1 = -1; double Z2K2L0 = 3.289868;
    double Z3K1L0 = 1.813376; double Z3K1L1 = -3.141593;
    double Z4K0L0 = 0.722126; double Z4K0L1 = -0.827216; double Z4K0L2 = 0.5; double Z4K2L0 = 0.696907; double Z4K2L1 = -3.289868; double Z4K4L0 = -2.164646;
    double Z5K1L0 = 2.268627; double Z5K1L1 = -2.598775; double Z5K1L2 = 1.570786; double Z5K3L0 = -3.776373;
    double Z6K0L0 = 0.730658; double Z6K0L1 = -0.991430; double Z6K0L2 = 0.538608; double Z6K0L3 = -0.166667; 
    double Z6K2L0 = 0.569598; double Z6K2L1 = -1.519374; double Z6K2L2 = 1.644934; double Z6K4L0 = -4.167149; double Z6K4L1 = 2.164646; double Z6K6L0 = 2.034686;

    double t1 = g.alpha*g.Z; double t2 = betaW/betaP; double t3=TMath::Log(betaP);
    double F2 = Z0K0L0*TMath::Power(t1,0)*TMath::Power(t2,0)*TMath::Power(t3,0)
            + Z1K1L0*TMath::Power(t1,1)*TMath::Power(t2,1)*TMath::Power(t3,0)
            + Z2K0L0*TMath::Power(t1,2)*TMath::Power(t2,0)*TMath::Power(t3,0)
            + Z2K0L1*TMath::Power(t1,2)*TMath::Power(t2,0)*TMath::Power(t3,1)
            + Z2K2L0*TMath::Power(t1,2)*TMath::Power(t2,2)*TMath::Power(t3,0)
            + Z3K1L0*TMath::Power(t1,3)*TMath::Power(t2,1)*TMath::Power(t3,0)
            + Z3K1L1*TMath::Power(t1,3)*TMath::Power(t2,1)*TMath::Power(t3,1)
            + Z4K0L0*TMath::Power(t1,4)*TMath::Power(t2,0)*TMath::Power(t3,0)
            + Z4K0L1*TMath::Power(t1,4)*TMath::Power(t2,0)*TMath::Power(t3,1)
            + Z4K0L2*TMath::Power(t1,4)*TMath::Power(t2,0)*TMath::Power(t3,2)
            + Z4K2L0*TMath::Power(t1,4)*TMath::Power(t2,2)*TMath::Power(t3,0)
            + Z4K2L1*TMath::Power(t1,4)*TMath::Power(t2,2)*TMath::Power(t3,1)
            + Z4K4L0*TMath::Power(t1,4)*TMath::Power(t2,4)*TMath::Power(t3,0)
            + Z5K1L0*TMath::Power(t1,5)*TMath::Power(t2,1)*TMath::Power(t3,0)
            + Z5K1L1*TMath::Power(t1,5)*TMath::Power(t2,1)*TMath::Power(t3,1)
            + Z5K1L2*TMath::Power(t1,5)*TMath::Power(t2,1)*TMath::Power(t3,2)
            + Z5K3L0*TMath::Power(t1,5)*TMath::Power(t2,3)*TMath::Power(t3,0)
            + Z6K0L0*TMath::Power(t1,6)*TMath::Power(t2,0)*TMath::Power(t3,0)
            + Z6K0L1*TMath::Power(t1,6)*TMath::Power(t2,0)*TMath::Power(t3,1)
            + Z6K0L2*TMath::Power(t1,6)*TMath::Power(t2,0)*TMath::Power(t3,2)
            + Z6K0L3*TMath::Power(t1,6)*TMath::Power(t2,0)*TMath::Power(t3,3)
            + Z6K2L0*TMath::Power(t1,6)*TMath::Power(t2,2)*TMath::Power(t3,0)
            + Z6K2L1*TMath::Power(t1,6)*TMath::Power(t2,2)*TMath::Power(t3,1)
            + Z6K2L2*TMath::Power(t1,6)*TMath::Power(t2,2)*TMath::Power(t3,2)
            + Z6K4L0*TMath::Power(t1,6)*TMath::Power(t2,4)*TMath::Power(t3,0)
            + Z6K4L1*TMath::Power(t1,6)*TMath::Power(t2,4)*TMath::Power(t3,1)
            + Z6K6L0*TMath::Power(t1,6)*TMath::Power(t2,6)*TMath::Power(t3,0);

    return F1*F2;
}


void fermiFunc::Plot() {
    TGraph* g0 = new TGraph();
    TGraph* g1 = new TGraph();
    const int nPoints = 1000;
    for(int iPoint=0; iPoint<nPoints; iPoint++) {
        g0->SetPoint(iPoint, g.E0/nPoints*(iPoint+1), phasespace(g.E0/nPoints*(iPoint+1)));
        g1->SetPoint(iPoint, g.E0/nPoints*(iPoint+1), phasespace(g.E0/nPoints*(iPoint+1))*fermiCorr(g.E0/nPoints*(iPoint+1)));
    }
    g0->SetLineColor(kGreen+1);
    g1->SetLineColor(kBlue);
    g0->SetName("phasespace");
    g1->SetName("fermi");


    TFile* graphFermi = new TFile("Fermi.root", "RECREATE");
    g0->Write();
    g1->Write();
    graphFermi->Close();
}


double gFermiFcn(double* x, double *p) {
    double betaE = x[0];
    if(gFermiFunc == 0 ) gFermiFunc = new fermiFunc();
    return gFermiFunc->fermiCorr(betaE)*gFermiFunc->phasespace(betaE);
}
