#include "junoChiFunction.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"
#include "junoParameters.hh"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

gammaResponse* junoChiFunction::gamma_array[9];
string junoChiFunction::gamma_name[9];
junoB12_simplified* junoChiFunction::b12data;

double junoChiFunction::m_chi2    = 0.;
double junoChiFunction::m_chi2Min = 100000;
double junoChiFunction::m_chi2B12 = 0;
double junoChiFunction::m_chi2Gam = 0;

int junoChiFunction::m_nParameter = 3;
double junoChiFunction::m_bestFit[20] = {0.};
double junoChiFunction::m_bestFitError[20] = {0.};

bool junoChiFunction::m_DoFit = false;
int junoChiFunction::m_nData;

bool junoChiFunction::m_doGamFit = junoParameters::fitGammaSources;
bool junoChiFunction::m_doB12Fit = junoParameters::fitB12;

junoChiFunction::junoChiFunction() {
    
    m_nData = 0;

    if (m_doGamFit) {
        gamma_array[m_nData] = new gammaResponse("Cs137", 100, 600, 1000);  gamma_name[m_nData] = "Cs137"; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Mn54", 100, 900, 1300);   gamma_name[m_nData] = "Mn54";  m_nData++;
        gamma_array[m_nData] = new gammaResponse("Ge68", 100, 1100, 1500);  gamma_name[m_nData] = "Ge68";  m_nData++;
        gamma_array[m_nData] = new gammaResponse("K40", 100, 1750, 2250);   gamma_name[m_nData] = "K40";   m_nData++;
        gamma_array[m_nData] = new gammaResponse("nH", 100, 2800, 3500);    gamma_name[m_nData] = "nH";    m_nData++;
        gamma_array[m_nData] = new gammaResponse("Co60", 100, 3200, 3700);  gamma_name[m_nData] = "Co60";  m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmBe", 100, 6000, 6800);  gamma_name[m_nData] = "AmBe";  m_nData++;
        gamma_array[m_nData] = new gammaResponse("nC12", 100, 6700, 7600);  gamma_name[m_nData] = "nC12";  m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmC", 100, 8400, 9400);   gamma_name[m_nData] = "AmC";   m_nData++;
    }

    if (m_doB12Fit) {
        b12data = new junoB12_simplified(100, 4500, 17500);
    }
}


junoChiFunction::~junoChiFunction()
{
    if (m_doGamFit) {
        for (int i=0; i<m_nData; i++) {
            delete gamma_array[i];
        }
    }

    if (m_doB12Fit) {
        delete b12data;
    }
}


void junoChiFunction::LoadData()
{
    if (m_doGamFit) {
        for(int i=0; i<m_nData; i++) {
            gamma_array[i]->LoadData();
        }
    }

    if (m_doB12Fit) {
        b12data->Initialize();
    }
}


double junoChiFunction::GetChi2( double maxChi2)
{
    double chi2 = 0;

    if (m_doGamFit) {
        for (int i=0; i<m_nData; i++) {
            chi2 += gamma_array[i]->GetChi2();
        }
    }

    if (m_doB12Fit) {
        chi2 += b12data->GetChi2();
    }

    cout << "current chi2 = " << chi2 << endl;

    return chi2;
}


void junoChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


void junoChiFunction::SetParameters(double *par)
{
    int iPar = 0;

    if (junoParameters::scintillatorParameterization == kSimulation) {
        electronQuench::setEnergyScale(par[iPar]);                  iPar++;
        electronQuench::setBirk1(par[iPar]);                        iPar++;
    }

    if (junoParameters::cerenkovMode == "kSimulationCer") {
        electronCerenkov::setkC(par[iPar]);                         iPar++;
    }

    if (junoParameters::pesigmaMode == "kTotal") {
        electronResponse::setra(par[iPar]);                         iPar++;
        electronResponse::setrb(par[iPar]);                         iPar++;
        electronResponse::setrc(par[iPar]);                         iPar++;
        electronResponse::SetParameters();
    }

    if (junoParameters::pesigmaMode == "kNPE") {
        electronResponse::setq0             (par[iPar]);            iPar++;
        electronResponse::setq1             (par[iPar]);            iPar++;
        electronResponse::setq2             (par[iPar]);            iPar++;
        electronResponse::SetParameters();

    }

    if (m_doGamFit) {
        for (int i=0; i<m_nData; i++) {
            gamma_array[i]->SetAmp(par[iPar]);                      iPar++;
        }
    }

}


double junoChiFunction::GetChiSquare(double maxChi2)
{
    junoMinuit = new TMinuit();
    junoMinuit->SetFCN(ChisqFCN);
    junoMinuit->SetPrintLevel(1);

    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    junoMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    if (junoParameters::scintillatorParameterization == kSimulation) {
        junoMinuit->mnparm(iPar, "scale", 1410, 1, 1300, 1500, ierrflag);   iPar++;
        junoMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-4, 5e-3, 7.5e-3, ierrflag);     iPar++;
    }


    if (junoParameters::cerenkovMode == "kSimulationCer") {
        junoMinuit->mnparm(iPar, "kC", 1.00, 0.001, 0, 1.5, ierrflag);            iPar++;
    }

    if (junoParameters::pesigmaMode == "kTotal") {
        junoMinuit->mnparm(iPar, "ra", 0, 0.01, -20, 20, ierrflag );              iPar++;
        junoMinuit->mnparm(iPar, "rb", 1315, 1, 1250, 1850, ierrflag);            iPar++;
        junoMinuit->mnparm(iPar, "rc", 160, 1, 100, 220, ierrflag);               iPar++;
    }

    if (junoParameters::pesigmaMode == "kNPE") {
        junoMinuit->mnparm(iPar, "q0", 22.333, 0.01, -100, 100, ierrflag);        iPar++;
        junoMinuit->mnparm(iPar, "q1", 0.983, 0.001, 0, 5, ierrflag);             iPar++;
        junoMinuit->mnparm(iPar, "q2", 6e-5, 1e-6, 1e-5, 1e-4, ierrflag);         iPar++;
    }

    if (m_doGamFit) {
        for(int j=0; j<m_nData; j++) {
            junoMinuit->mnparm(iPar, gamma_name[j].c_str(), 200, 1, 140, 300, ierrflag); iPar++;
        }
    }

    // Minimization strategy
    junoMinuit->SetErrorDef(1);
    arglist[0]=2;
    junoMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 5000; //maxCalls
    arglist[1] = 0.01; // tolerance
    junoMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    junoMinuit->fCstatu.Data();

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    junoMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    m_nParameter = junoMinuit->GetNumPars();
	for(int i=0; i<m_nParameter; i++)
	{
	    junoMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
	}

    m_DoFit = true;

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << " with nData = " << m_nData << " and nPar = " << m_nParameter << endl;
    cout << " ====================== " << endl;
    delete junoMinuit;
    return min;
}


void junoChiFunction::Plot()
{
    if (m_doGamFit) {
        for (int i=0; i<m_nData; i++) {
            gamma_array[i]->SaveHist();
        }
    }

    if (m_doB12Fit) {
        b12data->Plot();
    }
}

