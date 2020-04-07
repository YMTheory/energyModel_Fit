#include "gammaNLChiFunction.hh"
#include "gammaNLExperiment.hh"

#include <iostream>

using namespace std;

double gammaNLChiFunction::errGamma = 1.0;  // 100% error 

gammaNLExperiment* gammaNLChiFunction::mGammaNLExperiment = (gammaNLExperiment*)0;
gammaNLChiFunction::gammaNLChiFunction(gammaNLExperiment* aGammaNLExperiment)
{
    mGammaNLExperiment = aGammaNLExperiment;
}

gammaNLChiFunction::~gammaNLChiFunction()
{
    delete mGammaNLExperiment;
}

void gammaNLChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    double delta =0; double chisq = 0;

    cout << "kA " << par[0] << " kB " << par[1] << " kC " << par[2] << " err " << par[4]  << endl;

    // the first 3 parameters are model parameters ... 
    int varId = 3;

    double varGamma;
    varGamma = (par[varId]);

    mGammaNLExperiment->Calculate(par);
    TGraph* mFitgammaNL = mGammaNLExperiment->GetFitGammaNL();
    TGraphErrors* mTruegammaNL = mGammaNLExperiment->GetTrueGammaNL();

    const int nPoints = mTruegammaNL->GetN();
    double *Etrue = mTruegammaNL->GetX();
    double *nl_obs = mTruegammaNL->GetY();
    double *nl_err_obs = mTruegammaNL->GetEY();
    double *nl_pred = mFitgammaNL->GetY();  
    for(int iPoint=0; iPoint<nPoints; iPoint++){
        //delta = (nl_pred[iPoint] - nl_obs[iPoint]) / nl_err_obs[iPoint];
        cout << nl_pred[iPoint] << " " << nl_obs[iPoint] << " " << nl_err_obs[iPoint] << " " << varGamma << endl;
        delta = (nl_pred[iPoint]*(1+varGamma) - nl_obs[iPoint]) / nl_err_obs[iPoint];
        chisq += delta*delta;
    }

    chisq += ( varGamma / errGamma) * (varGamma/errGamma) ;
    cout << "chisq ===>  " << chisq << endl;
    fval = chisq;
}



double gammaNLChiFunction::GetChiSquare()
{
    gammaNLMinuit = new TMinuit();
    gammaNLMinuit->SetFCN(ChisqFCN);
    gammaNLMinuit->SetPrintLevel(1);
    
    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    gammaNLMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // Configurate parameters
    gammaNLMinuit->mnparm(iPar, "kA", 0.98, 0.001, 0.5, 1.5, ierrflag); iPar++;
    gammaNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 1e-3, 1e-2, ierrflag); iPar++;
    gammaNLMinuit->mnparm(iPar, "kC", 1.0, 0.001, 0.5, 1.5, ierrflag); iPar++;
    gammaNLMinuit->mnparm(iPar, "errGamma", 0, 0.1*errGamma, 0, 0, ierrflag); iPar++;
    
    //gammaNLMinuit->FixParameter(2);

    // Minimization strategy
    gammaNLMinuit->SetErrorDef(1);
    arglist[0]=2;
    gammaNLMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 100; //maxCalls
    arglist[1] = 0.01; // tolerance
    gammaNLMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    gammaNLMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    //gammaNLMinuit->SetErrorDef(1);
    //TGraph* graph1 = (TGraph*) gammaNLMinuit->Contour(50,0,1);
    //graph1->SetMarkerStyle(20);
    //graph1->SetMarkerColor(kBlue+1);
    //graph1->SetLineColor(kBlue+1);

    //TCanvas* cc = new TCanvas();
    //graph1->Draw("APL");
    //cc->SaveAs("tmp.png");

    delete gammaNLMinuit;
    return min;
}
