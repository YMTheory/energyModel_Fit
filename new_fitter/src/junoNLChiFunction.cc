#include "junoNLChiFunction.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"
#include "junoParameters.hh"
#include "junoSpectrum.hh"

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

gammaData* junoNLChiFunction::Cs137Data;
gammaData* junoNLChiFunction::Mn54Data;
gammaData* junoNLChiFunction::nHData;
gammaData* junoNLChiFunction::K40Data;
gammaData* junoNLChiFunction::Co60Data;
gammaData* junoNLChiFunction::Tl208Data;
gammaData* junoNLChiFunction::nC12Data;
gammaData* junoNLChiFunction::O16Data;
gammaData* junoNLChiFunction::nFe56Data;

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

int junoNLChiFunction::m_nGam;
int junoNLChiFunction::m_nData;
std::string junoNLChiFunction::source_name[20];
gammaData* junoNLChiFunction::gammaData_array[20];

bool junoNLChiFunction::m_doGamFit = junoParameters::fitGammaSources;
bool junoNLChiFunction::m_doB12Fit = junoParameters::fitB12;

string junoNLChiFunction::m_nonlMode;

junoNLChiFunction::junoNLChiFunction() {

    m_nData = 0;
    m_nGam  = 0;

    if (m_doGamFit) {

        Cs137Data = new gammaData("Cs137", 700, 1100, 100);
        source_name[m_nData] = "Cs137"; 
        gammaData_array[m_nData] = Cs137Data;
        m_nData++;
        m_nGam++;

        Mn54Data  = new gammaData("Mn54", 900, 1300, 100);
        source_name[m_nData] = "Mn54"; 
        gammaData_array[m_nData] = Mn54Data;
        m_nData++;
        m_nGam++;

        K40Data  = new gammaData("K40", 900, 1300, 100);
        source_name[m_nData] = "K40"; 
        gammaData_array[m_nData] = K40Data;
        m_nData++;
        m_nGam++;

        nHData  = new gammaData("nH", 900, 1300, 100);
        source_name[m_nData] = "nH"; 
        gammaData_array[m_nData] = nHData;
        m_nData++;
        m_nGam++;

        //Co60Data  = new gammaData("Co60", 900, 1300, 100);
        //source_name[m_nData] = "Co60"; 
        //gammaData_array[m_nData] = Co60Data;
        //m_nData++;

        Tl208Data  = new gammaData("Tl208", 900, 1300, 100);
        source_name[m_nData] = "Tl208"; 
        gammaData_array[m_nData] = Tl208Data;
        m_nData++;
        m_nGam++;

        nC12Data  = new gammaData("nC12", 900, 1300, 100);
        source_name[m_nData] = "nC12"; 
        gammaData_array[m_nData] = nC12Data;
        m_nData++;
        m_nGam++;

        O16Data  = new gammaData("O16", 900, 1300, 100);
        source_name[m_nData] = "O16"; 
        gammaData_array[m_nData] = O16Data;
        m_nData++;
        m_nGam++;

        nFe56Data  = new gammaData("nFe56", 900, 1300, 100);
        source_name[m_nData] = "nFe56"; 
        gammaData_array[m_nData] = nFe56Data;
        m_nData++;
        m_nGam++;

    }

    // Nonlinearity mode
    m_nonlMode = junoParameters::m_nonlMode;
    cout << "Nonlinearity formula form " << m_nonlMode << endl;

    if (m_doB12Fit) {
        junoB12 = new junoSpectrum(14400, 100, 3, 2,
                             0, 15, 0, 15, m_nonlMode, "B12");
    }

    electronResponse::FuncConstruct();
}

junoNLChiFunction::~junoNLChiFunction() {
    if (m_doGamFit) {
        delete Cs137Data;
        delete Mn54Data;
        delete nHData;
        delete K40Data;
        delete Co60Data;
        delete Tl208Data;
        delete nC12Data;
        delete O16Data;
        delete nFe56Data;
    }
    if(m_doB12Fit)
        delete junoB12;
}

void junoNLChiFunction::LoadData()
{
    if (m_doGamFit) {
        for (int i=0; i<m_nData; i++) {
            gammaData_array[i]->LoadData();
        }
    }
    if (m_doB12Fit) {
        junoB12->LoadData();
    }
}

double junoNLChiFunction::GetChi2( double maxChi2 )
{
    double chi2 = 0;

    if (m_doGamFit ) {
        for(int iSource=0; iSource<m_nData; iSource++) {
            chi2 += gammaData_array[iSource]->GetChi2();
        }
    }
    if(m_doB12Fit)
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
    if (junoParameters::scintillatorParameterization == kSimulation) { 
        if (m_nonlMode == "histogram") {
            electronQuench::setkA               (par[0]);
            electronQuench::setBirk1            (par[1]);
            electronCerenkov::setkC             (par[2]);
            electronCerenkov::setEnergyScale    (par[3]);
        }

        if (m_nonlMode == "analytic") {
            electronResponse::setp0(par[0]);
            electronResponse::setp1(par[1]);
            electronResponse::setp2(par[2]);
            electronResponse::setp3(par[3]);
            electronResponse::SetParameters();
        }
    }

    if (junoParameters::scintillatorParameterization == kSimulationCalc) {
        electronQuench::setkA               (par[0]);
        electronCerenkov::setkC             (par[1]);
        electronQuench::setEnergyScale      (par[2]);
        electronCerenkov::setEnergyScale    (par[2]);
    }
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
    if (junoParameters::scintillatorParameterization == kSimulation) {
        if (m_nonlMode == "histogram") {
            junoNLMinuit->mnparm(iPar, "kA", 1.00, 0.001, 0., 0., ierrflag); iPar++;
            junoNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.1e-3, 7.5e-3, ierrflag); iPar++;
            junoNLMinuit->mnparm(iPar, "kC", 1.0, 0.001, 0., 0., ierrflag); iPar++;
            junoNLMinuit->mnparm(iPar, "energyScale", 3300.371/2.223, 1, 0, 0, ierrflag); iPar++;
        }

        if (m_nonlMode == "analytic") {
            junoNLMinuit->mnparm(iPar, "p0", 1.025, 0.001, 0.9, 1.1, ierrflag);    iPar++;
            junoNLMinuit->mnparm(iPar, "p1", 0.1122,0.0001,  0, 0.2, ierrflag);    iPar++;
            junoNLMinuit->mnparm(iPar, "p2", 1.394, 0.001, 1.1, 4.0, ierrflag);    iPar++;
            junoNLMinuit->mnparm(iPar, "p3", 5.55e-4, 1e-5, 1e-5, 1e-2, ierrflag); iPar++;
        }
    }

    if (junoParameters::scintillatorParameterization == kSimulationCalc) {
        junoNLMinuit->mnparm(iPar, "kA", 1.00, 0.001, 0.0, 0.0, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "kC", 1.00, 0.001, 0.0, 0.0, ierrflag); iPar++;
        junoNLMinuit->mnparm(iPar, "energyScale", 3300.371/2.223, 1, 1400, 1600, ierrflag); iPar++;
    }

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

    junoNLMinuit->fCstatu.Data();

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    junoNLMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    m_nParameter = junoNLMinuit->GetNumPars();
	for(int i=0; i<m_nParameter; i++)
	{
	    junoNLMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
		//cout<<"curvalue: "<<curvalue<<"	curerror: "<<curerror<<endl;
	}

    m_DoFit = true;

    if (m_doB12Fit)
        m_nData += junoB12->getNData();

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << " with nData = " << m_nData << " and nPar = " << m_nParameter << endl;
    cout << " ====================== " << endl;
    delete junoNLMinuit;
    return min;
}



void junoNLChiFunction::Plot()
{
    //electronResponse::Plot();
    if(m_doGamFit)
        GammaPlot();
    if(m_doB12Fit)
        junoB12->Plot();
}


void junoNLChiFunction::GammaPlot()
{
    cout << " >>> Draw Gamma NL Fitting Results <<< " << endl;

    if (not m_DoFit) {
        cout << "Fitting has not been finished ...";
        return;
    }


    TGraphErrors* gNonlData = new TGraphErrors();
    TGraphErrors* gNonlCalc = new TGraphErrors();
    gNonlData->SetName("gNonlData");
    gNonlCalc->SetName("gNonlCalc");

    double fit_pars[3] = {m_bestFit[0], m_bestFit[1], m_bestFit[2]};
    SetParameters(fit_pars);
    int index = 0;
    for(int iData=0; iData<m_nGam; iData++) {
        std::string source = source_name[iData];
        gammaData* tmpGammaData = gammaData_array[iData];
        //tmpGammaData->calcGammaNPE();
        double tmp_E       = tmpGammaData->GetEtrue();
        double tmp_pred    = tmpGammaData->GetNonlPred();
        double tmp_data    = tmpGammaData->GetNonlData();
        double tmp_dataErr = tmpGammaData->GetNonlDataErr(); 
        //cout << tmp_E << " " << tmp_data << " " << tmp_pred << endl
        gNonlData->SetPoint(index, tmp_E, tmp_data);
        gNonlData->SetPointError(index, 0, tmp_dataErr);
        gNonlCalc->SetPoint(index, tmp_E, tmp_pred);

        index++;
    }

    gNonlData->SetMarkerStyle(20);
    gNonlData->SetMarkerColor(kBlue+1);
    gNonlData->SetLineColor(kBlue+1);
    gNonlData->SetLineWidth(2);
    gNonlData->SetMarkerSize(1.0);
    gNonlCalc->SetMarkerStyle(21);
    gNonlCalc->SetMarkerColor(kRed+1);
    gNonlCalc->SetMarkerSize(1.0);
    gNonlCalc->SetLineColor(kRed+1);
    gNonlCalc->SetLineWidth(2);

    double nom_pars[3] = {1.0, 1.0, 3300.371/2.223};
    SetParameters(nom_pars);
    TGraphErrors* gNom = new TGraphErrors();
    for(int iData=0; iData<m_nGam; iData++) {
        std::string source = source_name[iData];
        gammaData* tmpGammaData = gammaData_array[iData];
        tmpGammaData->calcGammaResponse();
        //tmpGammaData->calcGammaNPE();
        double tmp_E       = tmpGammaData->GetEtrue();
        double tmp_pred    = tmpGammaData->GetNonlPred();
        //cout << tmp_E << " " << tmp_data << " " << tmp_pred << endl
        gNom->SetPoint(index, tmp_E, tmp_pred);

        index++;
    }
    gNom->SetMarkerStyle(21);
    gNom->SetMarkerColor(kOrange+1);
    gNom->SetMarkerSize(1.0);
    gNom->SetLineColor(kOrange+1);
    gNom->SetLineWidth(2);


    TCanvas* c1 = new TCanvas("Nonlinearity", "Nonlinearity");
    c1->cd(); c1->SetGrid();
    gNonlData->SetTitle("Nonlinearity Fitting; Etrue/MeV; Nonlinearity");
    //gNonlData->GetYaxis()->SetRangeUser(0.01,0.045);
    gNonlData->Draw("APL");
    gNonlCalc->Draw("P SAME");
    gNom->Draw("P SAME");
    TLegend* led = new TLegend();
    led->SetFillColor(kWhite);
    led->AddEntry(gNonlData, "data", "PL");
    led->AddEntry(gNonlCalc, "calc", "PL");
    led->Draw("SAME");

    c1->SaveAs("GamNLFit.root");    
    
}



