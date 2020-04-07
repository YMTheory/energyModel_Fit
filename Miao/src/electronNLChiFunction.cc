#include "electronNLChiFunction.hh"
#include "electronNLExperiment.hh"

#include <TMinuit.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>

#include <vector>
#include <iostream>

using namespace std;
double electronNLChiFunction::errElectron = 1.0;  // 100% error 

electronNLExperiment* electronNLChiFunction::melectronNLExperiment = (electronNLExperiment*)0;

electronNLChiFunction::electronNLChiFunction(electronNLExperiment* aelectronNLExperiment)
{
    melectronNLExperiment = aelectronNLExperiment;
}

electronNLChiFunction::~electronNLChiFunction()
{;}

void electronNLChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    double delta =0; double chisq = 0;

    // the first 3 parameters are model parameters ... 
    int varId = 3;

    double varElectron;
    varElectron = (par[varId]);

    melectronNLExperiment->Calculate(par);
    TGraph* mFitElectronNL = melectronNLExperiment->GetFitElectronNL();
    TGraphErrors* mTrueElectronNL = melectronNLExperiment->GetTrueElectronNL();

    const int nPoints = mTrueElectronNL->GetN();
    double *Etrue = mTrueElectronNL->GetX();
    double *nl_obs = mTrueElectronNL->GetY();
    double *nl_err_obs = mTrueElectronNL->GetEY();
    double *nl_pred = mFitElectronNL->GetY();  
    for(int iPoint=0; iPoint<nPoints; iPoint++){
        delta = (nl_pred[iPoint] - nl_obs[iPoint]) / nl_err_obs[iPoint];
        //delta = (nl_pred[iPoint] - nl_obs[iPoint]*(1+varElectron)) / nl_err_obs[iPoint];
        chisq += delta*delta;
    }

    //chisq += ( varElectron / errElectron) * (varElectron/errElectron) ;
    fval = chisq;
}

double electronNLChiFunction::GetChiSquare()
{
    electronNLMinuit = new TMinuit();
    electronNLMinuit->SetFCN(ChisqFCN);
    electronNLMinuit->SetPrintLevel(1);
    
    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    electronNLMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // Configurate parameters
    electronNLMinuit->mnparm(iPar, "kA", 0.98, 0.01, 0, 1, ierrflag); iPar++;
    electronNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 1e-3, 1e-2, ierrflag); iPar++;
    electronNLMinuit->mnparm(iPar, "kC", 1.0, 0.01, 0, 1, ierrflag); iPar++;
    //electronNLMinuit->mnparm(iPar, "errElectron", 0.0,  0.01, 0, 0, ierrflag); iPar++;
    
    //electronNLMinuit->FixParameter(2);

    // Minimization strategy
    electronNLMinuit->SetErrorDef(1);
    arglist[0]=2;
    electronNLMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 5000; //maxCalls
    arglist[1] = 0.01; // tolerance
    electronNLMinuit->mnexcm("MIGrad", arglist, 2, ierrflag);

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    electronNLMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    //electronNLMinuit->SetErrorDef(1);
    //TGraph* graph1 = (TGraph*) electronNLMinuit->Contour(10,0,1);
    //graph1->SetMarkerStyle(20);
    //graph1->SetMarkerColor(kBlue+1);
    //graph1->SetLineColor(kBlue+1);

    //TCanvas* cc = new TCanvas();
    //graph1->Draw("APL");
    //cc->SaveAs("tmp.png");

    delete electronNLMinuit;
    return min;
}

