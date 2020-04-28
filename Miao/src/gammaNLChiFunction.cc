#include "gammaNLChiFunction.hh"
#include "gammaNLExperiment.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

#include <iostream>

using namespace std;

double gammaNLChiFunction::errGamma = 1.0;  // 100% error 
double gammaNLChiFunction::m_chi2 = 0.;

int gammaNLChiFunction::m_nParameter = 4;
double gammaNLChiFunction::m_bestFit[20] = {0.};
double gammaNLChiFunction::m_bestFitError[20] = {0.};

bool gammaNLChiFunction::m_DoFit = false;

gammaNLChiFunction::gammaNLChiFunction()
{;}

gammaNLChiFunction::~gammaNLChiFunction()
{;}

double gammaNLChiFunction::GetChi2 ( double maxChi2 ) {
    m_chi2 = 0;
    m_chi2 += gammaNLExperiment::GetChi2(0);
    m_chi2 += (gammaNLExperiment::m_gammaScale/junoParameters::m_gammaError, 2);
    if(maxChi2>0 and m_chi2>maxChi2) return maxChi2;

    return m_chi2;
}


void gammaNLChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


double gammaNLChiFunction::GetChiSquare(double maxChi2)
{
    gammaNLMinuit = new TMinuit();
    gammaNLMinuit->SetFCN(ChisqFCN);
    gammaNLMinuit->SetPrintLevel(1);
    
    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    gammaNLMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // Configurate parameters
    gammaNLMinuit->mnparm(iPar, "kA", 0.98, 0.01, 0., 2.0, ierrflag); iPar++;
    gammaNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-4, 1e-4, 1e-2, ierrflag); iPar++;
    gammaNLMinuit->mnparm(iPar, "kC", 1.0, 0.01, 0.0, 2.0, ierrflag); iPar++;
    gammaNLMinuit->mnparm(iPar, "errGamma", 0, 0.1*errGamma, 0, 0, ierrflag); iPar++;
    
    //gammaNLMinuit->FixParameter(1);

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

	for(int i=0; i<m_nParameter; i++)
	{
	    gammaNLMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
		//cout<<"curvalue: "<<curvalue<<"	curerror: "<<curerror<<endl;
	}

    m_DoFit = true;

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;
    delete gammaNLMinuit;
    return min;
}


void gammaNLChiFunction::SetParameters(double *par) {
    electronQuench::setkA             (par[0]);
    electronQuench::setBirk1          (par[1]);
    electronCerenkov::setkC           (par[2]);
    gammaNLExperiment::setGammaScale  (par[3]);
}


void gammaNLChiFunction::Plot ()
{
    if(!m_DoFit) {
        cout << " >>> Have Not Done Fitting Yet ! <<< " << endl; return;
    }
    int nPoints = m_nParameter;
    double par[nPoints];
    for(int iPoint=0; iPoint<nPoints; iPoint++) {
        par[iPoint] = m_bestFit[iPoint];
    }
    SetParameters(par);
    gammaNLExperiment::Plot();
}
