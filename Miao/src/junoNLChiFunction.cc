#include "junoNLChiFunction.hh"
#include "electronNLExperiment.hh"
#include "gammaNLExperiment.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "BetaPrediction.hh"
#include "junoB12Data.hh"
#include "junoC11Data.hh"
#include "junoParameters.hh"

#include <iostream>

using namespace std;

junoB12Data* junoNLChiFunction::m_b12Data = 0;
junoC11Data* junoNLChiFunction::m_c11Data = 0;

double junoNLChiFunction::m_chi2    = 0.;
double junoNLChiFunction::m_chi2Min = 100000;
double junoNLChiFunction::m_chi2B12 = 0;
double junoNLChiFunction::m_chi2C11 = 0;

int junoNLChiFunction::m_nParameter = 3;
double junoNLChiFunction::m_bestFit[20] = {0.};
double junoNLChiFunction::m_bestFitError[20] = {0.};

bool junoNLChiFunction::m_DoFit = false;

junoNLChiFunction::junoNLChiFunction() {
    m_b12Data = new junoB12Data();
    m_b12Data->SetParameters();
    m_c11Data = new junoC11Data();
    m_c11Data->SetParameters();
}

junoNLChiFunction::~junoNLChiFunction() {
    delete m_b12Data;
}

void junoNLChiFunction::LoadData()
{
	std::cout << " ---> Start loading data " << std::endl;
    m_b12Data   ->LoadData(junoParameters::B12DataFile);
    m_c11Data   ->LoadData("./data/electron/C11.root");
}


double junoNLChiFunction::GetChi2(double maxChi2) {
    m_chi2 = 0;
    m_chi2B12 = 0;
    m_chi2C11 = 0;
    //m_chi2 += electronNLExperiment::GetChi2(0);

    if(junoParameters::fitB12) {
        m_chi2B12 = m_b12Data ->GetChi2();
        m_chi2 += m_chi2B12;
    }

    if(junoParameters::fitC11) {
        m_chi2C11 = m_c11Data ->GetChi2();
        m_chi2 += m_chi2C11;
    }
    m_chi2 += gammaNLExperiment::GetChi2(0);

    if(maxChi2>0 and m_chi2>maxChi2) return maxChi2;
    return m_chi2;
}


void junoNLChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


void junoNLChiFunction::SetParameters(double* par)
{
    //BetaPrediction::setK                (par[0]);   // B12 Spec Normalization
    electronQuench::setkA               (par[0]);
    electronQuench::setBirk1            (par[1]);
    electronCerenkov::setkC             (par[2]);
    gammaNLExperiment::setGammaScale    (par[3]);
    //electronQuench::setkA               ((1-par[2]*58.517/1481.06)/0.9796);
}



double junoNLChiFunction::GetChiSquare(double maxChi2)
{
    junoNLMinuit = new TMinuit();
    junoNLMinuit->SetFCN(ChisqFCN);
    junoNLMinuit->SetPrintLevel(1);
    
    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    junoNLMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // Configurate parameters
    junoNLMinuit->mnparm(iPar, "kA", 0.98, 0.01, 0.7, 1.2, ierrflag); iPar++;
    //junoNLMinuit->mnparm(iPar, "K", 0.116, 0.1,1.13, 0.001, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-4, 5e-3, 7e-3, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "kC", 1.0, 0.01, 0.7, 1.2, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "errGamma", 0, 0.01, 0, 0, ierrflag); iPar++;
    
    //junoNLMinuit->FixParameter(0);
    //junoNLMinuit->FixParameter(1);
    //junoNLMinuit->FixParameter(2);
    //junoNLMinuit->FixParameter(3);

    // Minimization strategy
    junoNLMinuit->SetErrorDef(1);
    arglist[0]=2;
    junoNLMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 1000; //maxCalls
    arglist[1] = 0.01; // tolerance
    junoNLMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    junoNLMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

	for(int i=0; i<m_nParameter; i++)
	{
	    junoNLMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
		//cout<<"curvalue: "<<curvalue<<"	curerror: "<<curerror<<endl;
	}

    m_DoFit = true;

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;
    delete junoNLMinuit;
    return min;
}


void junoNLChiFunction::Plot()
{
    if(!m_DoFit) {
        cout << " >>> Fitting has not been finished !<<< " << endl; return;
    }
    
    int nPoints = m_nParameter;
    double par[nPoints];
    for(int iPoint=0; iPoint<nPoints; iPoint++) {
        par[iPoint] = m_bestFit[iPoint];
    }
    SetParameters(par);
    //gammaNLExperiment::Plot();
    if(junoParameters::fitB12) {
        m_b12Data->Plot();
    }
    if(junoParameters::fitC11) {
        m_c11Data->Plot(); // 
    }
}







