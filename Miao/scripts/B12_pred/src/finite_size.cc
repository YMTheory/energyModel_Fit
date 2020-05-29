#include <global.hh>
#include <finite_size.hh>

#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>

finite_size* gFiniteSize = (finite_size*)0;

finite_size::finite_size()
{
    mFiniteSize = new TF1("finitesize_fcn", "", 0, g.E0);
}

double finite_size::elecmagCorr(double betaE)
{
    double betaW = betaE/0.511+1;
    double t1 = g.alpha*g.Z;
    double aminus1 = 0.115*t1 - 1.8123*TMath::Power(t1,2) +8.2498*TMath::Power(t1,3)-11.223*TMath::Power(t1,4)-14.854*TMath::Power(t1,5)+32.086*TMath::Power(t1,6);
    
    //electromagnetic finite_size corrections:
    double L0 = 1 + 13*(g.alpha*g.R)*(g.alpha*g.R)/60. - betaW*g.R*g.alpha*g.Z*(41-26*g.gamma)/(15*(2*g.gamma-1)) - g.alpha*g.Z*g.R*g.gamma*(17-2*g.gamma)/(30*betaW*(2*g.gamma-1)) + aminus1*g.R/betaW;
    
    return L0;
}

double finite_size::weakCorr(double betaE)
{
    double betaW = betaE/0.511+1;
    double t1 = g.alpha*g.Z;
    // weak-interaction finite-size correction: 
    double C0 = -233/630.*t1*t1 - (g.W0*g.R)*(g.W0*g.R)/5 + 2./35*g.W0*g.R*g.alpha*g.Z;
    double C1 = -21./35*g.R*t1 + 4./9*g.W0*g.R*g.R;
    double C2 = -4./9*g.R*g.R;
    double C = 1 + C0 + C1*betaW + C2*betaW*betaW;
    
    return C;
}

void finite_size::Plot()
{
    //TCanvas* cc = new TCanvas(); cc->cd();

    TGraph* graphElecmag = new TGraph();
    TGraph* graphWeak = new TGraph();
    const int nPoints = 1000;
    for(int iPoint=0; iPoint<nPoints; iPoint++) {
        double E = g.E0/nPoints*(iPoint+1);
        graphElecmag->SetPoint(iPoint, E, elecmagCorr(E) );
        graphWeak->SetPoint(iPoint, E, weakCorr(E));
    }
    graphElecmag->SetName("ElecFiniteSize");
    graphElecmag->SetLineColor(kBlue+1);
    graphWeak->SetName("WeakFiniteSize");
    graphWeak->SetLineColor(kGreen+1);
    //graphElecmag->Draw("AL");
    //graphWeak->Draw("L SAME");
        
    //TLegend* ll = new TLegend();
    //ll->AddEntry(graphElecmag, "electromagnetic finite-size", "l");
    //ll->AddEntry(graphWeak, "weak finite_size", "l");
    //ll->Draw("SAME");

    TFile* file = new TFile("FiniteSize.root", "RECREATE");
    graphElecmag->Write();
    graphWeak->Write();
    file->Close();
}

double gFiniteSizeFcn(double *x, double *p) {
    double betaE = x[0];
    if(gFiniteSize == 0) gFiniteSize = new finite_size();
    return gFiniteSize->elecmagCorr(betaE)*gFiniteSize->weakCorr(betaE);
}
