#include "electronResChiFunction.hh"
#include "electronResExperiment.hh"
#include "electronResol.hh"

#include "TMinuit.h"

#include <iostream>

using namespace std;

double electronResChiFunction::m_chi2 = 0.;
double electronResChiFunction::m_chi2Min = 10000;

int electronResChiFunction::m_nParameter = 3;
double electronResChiFunction::m_bestFit[20] = {0.};
double electronResChiFunction::m_bestFitError[20] = {0.};

bool electronResChiFunction::m_DoFit = false;

electronResChiFunction::electronResChiFunction()
{;}

electronResChiFunction::~electronResChiFunction()
{;}

double electronResChiFunction::GetChi2(double maxChi2) {
    m_chi2 = 0;
    m_chi2 += electronResExperiment::GetChi2(); 
    if(maxChi2>0 and m_chi2>maxChi2) return maxChi2;

    return m_chi2;
}

void electronResChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2(); 
}


double electronResChiFunction::GetChiSquare( double maxChi2 )
{
    electronResMinuit = new TMinuit();
    electronResMinuit->SetFCN(ChisqFCN);
    electronResMinuit->SetPrintLevel(1);
    
    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    electronResMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // Configurate parameters
    electronResMinuit->mnparm(iPar, "pA", 0.01, 0.001, 0, 0.1, ierrflag);           iPar++;
    electronResMinuit->mnparm(iPar, "pB", 0.0068, 0.0001, 0., 0.01, ierrflag);   iPar++;
    electronResMinuit->mnparm(iPar, "pC", 0.0, 0.001, 0, 1, ierrflag);            iPar++;
    //electronResMinuit->mnparm(iPar, "errElectron", 0.0,  0.01, 0, 0, ierrflag); iPar++;
    
    //electronResMinuit->FixParameter(0);
    //electronResMinuit->FixParameter(1);
    //electronResMinuit->FixParameter(2);

	m_nParameter = electronResMinuit->GetNumPars();

    // Minimization strategy
    electronResMinuit->SetErrorDef(1);
    arglist[0]=2;
    electronResMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 50000; //maxCalls
    arglist[1] = 0.01; // tolerance
    electronResMinuit->mnexcm("MIGrad", arglist, 2, ierrflag);

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    electronResMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

	for(int i=0; i<m_nParameter; i++)
	{
	    electronResMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
		//cout<<"curvalue: "<<curvalue<<"	curerror: "<<curerror<<endl;
	}


    m_DoFit = true;

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;
    delete electronResMinuit;
    m_chi2Min = min;
    return min;
}


void electronResChiFunction::SetParameters(double *par)
{
    electronResol::setpA(par[0]);
    electronResol::setpB(par[1]);
    electronResol::setpC(par[2]);
}


void electronResChiFunction::Plot() {
    if(!m_DoFit) {
        cout << " >>> Have Not Done Fitting Yet ! <<< " << endl; return;
    } else {
        int nPars = m_nParameter;
        double par[nPars];
        for(int iPar = 0; iPar<nPars; iPar++) {
            par[iPar] = m_bestFit[iPar];
        }
        SetParameters(par);
        electronResExperiment::PlotB12();
        electronResExperiment::PlotC11();
    }
}

