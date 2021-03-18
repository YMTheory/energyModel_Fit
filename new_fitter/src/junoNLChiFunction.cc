#include "junoNLChiFunction.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "junoParameters.hh"
#include "junoSpectrum.hh"

using namespace std;

junoSpectrum* junoNLChiFunction::junoB12;

double junoNLChiFunction::m_chi2    = 0.;
double junoNLChiFunction::m_chi2Min = 100000;
double junoNLChiFunction::m_chi2B12 = 0;
double junoNLChiFunction::m_chi2C11 = 0;
double junoNLChiFunction::m_chi2C10 = 0;
double junoNLChiFunction::m_chi2Gam = 0;

int junoNLChiFunction::m_nParameter = 3;
double junoNLChiFunction::m_bestFit[20] = {0.};
double junoNLChiFunction::m_bestFitError[20] = {0.};

bool junoNLChiFunction::m_DoFit = false;

double junoNLChiFunction::final_kA = 0;
double junoNLChiFunction::final_kB = 0;
double junoNLChiFunction::final_kC = 0;

junoNLChiFunction::junoNLChiFunction() {
    junoB12 = new junoSpectrum(14400, 100, 3, 1,
                               0, 15, 0, 15, "B12");
}

junoNLChiFunction::~junoNLChiFunction() {
    delete junoB12;
}

void junoNLChiFunction::LoadData()
{
    junoB12->LoadData();
}

double junoNLChiFunction::GetChi2( double maxChi2 )
{
    double chi2 = 0;
    chi2 += junoB12->GetChi2();

    return chi2;
}

void junoNLChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}

void junoNLChiFunction::SetParameters(double *par)
{
    electronQuench::setkA               (par[0]);
    electronQuench::setBirk1            (par[1]);
    electronCerenkov::setkC             (par[2]);
    electronCerenkov::setEnergyScale    (par[3]);
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
    junoNLMinuit->mnparm(iPar, "kA", 0.96, 0.001, 0.9, 1.0, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.5e-3, 7.5e-3, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "kC", 1.0, 0.001, 0.9, 1.20, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "energyScale", 1500, 1, 1400, 1600, ierrflag); iPar++;
    
    //junoNLMinuit->FixParameter(0);
    //junoNLMinuit->FixParameter(1);
    //junoNLMinuit->FixParameter(2);
    //junoNLMinuit->FixParameter(3);

    // Minimization strategy
    junoNLMinuit->SetErrorDef(1);
    arglist[0]=2;
    junoNLMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 5000; //maxCalls
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
    junoB12->Plot();
}





