#include "BetaPrediction.hh"
#include "junoParameters.hh"
#include <TMath.h>
#include <TGraph.h>
#include <TFile.h>

using namespace std;

double BetaPrediction::gamma = TMath::Sqrt(1-(alpha*Z)*(alpha*Z));
double BetaPrediction::R     = 0.0029*TMath::Power(A,1/3)+0.0063*TMath::Power(A,-1/3)-0.017*TMath::Power(A,-1);

double BetaPrediction::predSpec(double betaE) {
    double betaW = betaE/0.511+1;
    double betaP = TMath::Sqrt(betaW*betaW-1);
    double sp = betaP*betaP*(W0-betaW)*(W0-betaW);
    double F1 = 2*(gamma+1)*TMath::Power(2*R, 2*(gamma-1))*0.253564; // energy-independent part
    // Table for expansions ...
    double Z0K0L0 = 1; 
    double Z1K1L0 = 3.141593;
    double Z2K0L0 = 0.577216; double Z2K0L1 = -1; double Z2K2L0 = 3.289868;
    double Z3K1L0 = 1.813376; double Z3K1L1 = -3.141593;
    double Z4K0L0 = 0.722126; double Z4K0L1 = -0.827216; double Z4K0L2 = 0.5; double Z4K2L0 = 0.696907; double Z4K2L1 = -3.289868; double Z4K4L0 = -2.164646;
    double Z5K1L0 = 2.268627; double Z5K1L1 = -2.598775; double Z5K1L2 = 1.570786; double Z5K3L0 = -3.776373;
    double Z6K0L0 = 0.730658; double Z6K0L1 = -0.991430; double Z6K0L2 = 0.538608; double Z6K0L3 = -0.166667; 
    double Z6K2L0 = 0.569598; double Z6K2L1 = -1.519374; double Z6K2L2 = 1.644934; double Z6K4L0 = -4.167149; double Z6K4L1 = 2.164646; double Z6K6L0 = 2.034686;

    double t1 = alpha*Z; double t2 = betaW/betaP; double t3=TMath::Log(betaP);
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

    double aminus1 = 0.115*t1 - 1.8123*TMath::Power(t1,2) +8.2498*TMath::Power(t1,3)-11.223*TMath::Power(t1,4)-14.854*TMath::Power(t1,5)+32.086*TMath::Power(t1,6);
    
    //electromagnetic finite_size corrections:
    double L0 = 1 + 13*(alpha*R)*(alpha*R)/60. - betaW*R*alpha*Z*(41-26*gamma)/(15*(2*gamma-1)) - alpha*Z*R*gamma*(17-2*gamma)/(30*betaW*(2*gamma-1)) + aminus1*R/betaW;
    // weak-interaction finite-size correction: 
    double C0 = -233/630.*t1*t1 - (W0*R)*(W0*R)/5 + 2./35*W0*R*alpha*Z;
    double C1 = -21./35*R*t1 + 4./9*W0*R*R;
    double C2 = -4./9*R*R;
    double C = 1 + C0 + C1*betaW + C2*betaW*betaW;

    //screening effect
    double Z_bar = Z;
    double V0 = alpha*alpha*TMath::Power(Z_bar, 0.75)*1.4584000;  // linear interpolation from table
    double betaW_bar = betaW - V0;
    double betaP_bar = TMath::Sqrt(betaW_bar*betaW_bar-1);
    double y = alpha*Z*betaW/betaP;
    double y_bar = alpha*Z*betaW_bar/betaP_bar;
    double S1 = betaW_bar/betaW * TMath::Power(betaP_bar/betaP, 2*gamma-1)*TMath::Exp(TMath::Pi()*(y_bar-y)) * 0.253564;
    double Z0K0 = 1;
    double Z2K0 = 0.577216; double Z2K2 = -1.644934;
    double Z4K0 = 0.722126; double Z4K2 = -2.151539; double Z4K4 = 1.894066;
    double Z6K0 = 0.730658; double Z6K2 = -2.993953; double Z6K4 = 4.107516; double Z6K6 = -1.971102;
    double tt1 = alpha*Z; double tt2 = betaW_bar/betaP_bar;
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

    // weak_magnetism correction
    double delta = 0.0048;

    return K*sp*F1*F2*L0*C*S*(1+delta*betaW);

}

void BetaPrediction::Plot() {
    TGraph* graphPred = new TGraph();
    for(int ibin=0; ibin<m_bin; ibin++ ) {
        double E = E0/m_bin*(ibin+1);
        graphPred->SetPoint(ibin, E, predSpec(E));
    }

    graphPred->SetLineColor(kBlue);
    graphPred->SetMarkerColor(kBlue);
    graphPred->SetMarkerSize(0.3);
    graphPred->SetName("B12Calc");

    TFile* file = new TFile(junoParameters::B12CalcFile.c_str(), "recreate");
    graphPred->Write();
    file->Close();
}

