#include "junoNLChiFunction.hh"
#include "gammaResol.hh"
#include "gammaNLExperiment.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "junoB12Data.hh"
#include "junoC11Data.hh"
#include "junoC10Data.hh"
#include "junoParameters.hh"

#include <iostream>
#include <fstream>

#include "TGraphErrors.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

gammaResol* junoNLChiFunction::Cs137Data = 0;
gammaResol* junoNLChiFunction::Mn54Data = 0;
gammaResol* junoNLChiFunction::K40Data = 0;
gammaResol* junoNLChiFunction::nHData = 0;
gammaResol* junoNLChiFunction::Co60Data = 0;
gammaResol* junoNLChiFunction::Tl208Data = 0;
gammaResol* junoNLChiFunction::nC12Data = 0;
gammaResol* junoNLChiFunction::O16Data = 0;
gammaResol* junoNLChiFunction::nFe56Data = 0;

junoB12Data* junoNLChiFunction::m_b12Data = 0;
junoC11Data* junoNLChiFunction::m_c11Data = 0;
junoC10Data* junoNLChiFunction::m_c10Data = 0;

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
bool junoNLChiFunction::m_gridSearch = false;

double junoNLChiFunction::final_kA = 0;
double junoNLChiFunction::final_kB = 0;
double junoNLChiFunction::final_kC = 0;
double junoNLChiFunction::m_kALimit = 0.01;
double junoNLChiFunction::m_kBLimit = 3e-4;
double junoNLChiFunction::m_kCLimit = 0.01;

int junoNLChiFunction::m_nData = 0;
std::map<std::string, gammaResol*> junoNLChiFunction::mapGammaResol;
std::string junoNLChiFunction::source_name[20];

junoNLChiFunction::junoNLChiFunction() {

    Cs137Data = new gammaResol("Cs137", 700, 1100, 100);
    Cs137Data->LoadData(); 
    source_name[m_nData] = "Cs137"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Cs137", Cs137Data));

    Mn54Data  = new gammaResol("Mn54", 900, 1300, 100);
    Mn54Data->LoadData(); 
    source_name[m_nData] = "Mn54"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Mn54", Mn54Data));

    K40Data  = new gammaResol("K40", 900, 1300, 100);
    K40Data->LoadData(); 
    source_name[m_nData] = "K40"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("K40", K40Data));

    nHData  = new gammaResol("nH", 900, 1300, 100);
    nHData->LoadData(); 
    source_name[m_nData] = "nH"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("nH", nHData));

    Co60Data  = new gammaResol("Co60", 900, 1300, 100);
    Co60Data->LoadData(); 
    source_name[m_nData] = "Co60"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Co60", Co60Data));

    Tl208Data  = new gammaResol("Tl208", 900, 1300, 100);
    Tl208Data->LoadData(); 
    source_name[m_nData] = "Tl208"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Tl208", Tl208Data));

    nC12Data  = new gammaResol("nC12", 900, 1300, 100);
    nC12Data->LoadData(); 
    source_name[m_nData] = "nC12"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("nC12", nC12Data));

    O16Data  = new gammaResol("O16", 900, 1300, 100);
    O16Data->LoadData(); 
    source_name[m_nData] = "O16"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("O16", O16Data));

    nFe56Data  = new gammaResol("nFe56", 900, 1300, 100);
    nFe56Data->LoadData(); 
    source_name[m_nData] = "nFe56"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("nFe56", nFe56Data));

    // continuous beta spectrum
    m_b12Data = new junoB12Data();
    m_b12Data->SetParameters();
    m_c11Data = new junoC11Data();
    m_c11Data->SetParameters();
    m_c10Data = new junoC10Data();
    m_c10Data->SetParameters();
}

junoNLChiFunction::~junoNLChiFunction() {
    delete m_b12Data;
    delete m_c11Data;
    delete m_c10Data;
}

void junoNLChiFunction::LoadData()
{
	std::cout << " ---> Start loading data " << std::endl;
    m_b12Data   ->LoadData(junoParameters::B12DataFile);
    m_c11Data   ->LoadData(junoParameters::C11DataFile);
    m_c10Data   ->LoadData(junoParameters::C10DataFile);
}


double junoNLChiFunction::GetChi2(double maxChi2) {
    m_chi2    = 0;
    m_chi2Gam = 0;
    m_chi2B12 = 0;
    m_chi2C11 = 0;
    m_chi2C10 = 0;
    //m_chi2 += electronNLExperiment::GetChi2(0);

    if(junoParameters::fitGammaSources) {
        m_chi2Gam += Cs137Data->GetChi2();
        m_chi2Gam += Mn54Data->GetChi2();
        m_chi2Gam += K40Data->GetChi2();
        m_chi2Gam += nHData->GetChi2();
        m_chi2Gam += Co60Data->GetChi2();
        m_chi2Gam += Tl208Data->GetChi2();
        m_chi2Gam += nC12Data->GetChi2();
        m_chi2Gam += O16Data->GetChi2();
        m_chi2Gam += nFe56Data->GetChi2();

        m_chi2 += m_chi2Gam;
    }

    if(junoParameters::fitB12) {
        m_chi2B12 = m_b12Data ->GetChi2();
        m_chi2 += m_chi2B12;
    }

    if(junoParameters::fitC11) {
        m_chi2C11 = m_c11Data ->GetChi2();
        m_chi2 += m_chi2C11;
    }

    if(junoParameters::fitC10) {
        m_chi2C10 = m_c10Data ->GetChi2();
        m_chi2 += m_chi2C10;
    }

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
    junoSpectrum::setGammaScale         (par[3]);
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
    junoNLMinuit->mnparm(iPar, "kA", 0.96, 0.001, 0.9, 1.0, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.5e-3, 7.5e-3, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "kC", 1.0, 0.001, 0.9, 1.00, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "errGamma", 0, 0.01, 0, 0, ierrflag); iPar++;
    
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
    if(!m_DoFit) {
        cout << " >>> Fitting has not been finished !<<< " << endl; return;
    }
    
    if(m_gridSearch) {
        int nPoints = 3;
        double par[nPoints];
        par[0] = final_kA; par[1] = final_kB; par[2] = final_kC;
        SetParameters(par);
    } else {
        int nPoints = m_nParameter;
        double par[nPoints];
        for(int iPoint=0; iPoint<nPoints; iPoint++) {
            par[iPoint] = m_bestFit[iPoint];
        }
        SetParameters(par);
    }

    if(junoParameters::fitGammaSources) {
        GammaPlot();
    }
    if(junoParameters::fitB12) {
        m_b12Data->Plot();
    }
    if(junoParameters::fitC11) {
        m_c11Data->Plot(); // 
    }

    if(junoParameters::fitC10) {
        m_c10Data->Plot(); // 
    }

}


void junoNLChiFunction::GammaPlot() {
    TGraphErrors* gNonlData = new TGraphErrors();
    TGraph* gNonlCalc = new TGraph();
    gNonlData->SetName("gNonlData");
    gNonlCalc->SetName("gNonlCalc");
    int index = 0;
    for(int iData=0; iData<m_nData; iData++) {
        std::string source = source_name[iData];
        gammaResol* tmpGammaData = mapGammaResol.find(source)->second;
        //tmpGammaData->calcGammaNPE();
        double tmp_E       = tmpGammaData->GetEtrue();
        double tmp_pred    = tmpGammaData->GetNonlPred();
        double tmp_data    = tmpGammaData->GetNonlData();
        double tmp_dataErr = tmpGammaData->GetNonlDataErr(); 
        //cout << tmp_E << " " << tmp_data << " " << tmp_pred << endl;
        gNonlData->SetPoint(index, tmp_E, tmp_data);
        gNonlData->SetPointError(index, 0, tmp_dataErr);
        gNonlCalc->SetPoint(index, tmp_E, tmp_pred);
        index++;
    }

    gNonlData->SetMarkerStyle(20);
    gNonlData->SetMarkerColor(kBlue+1);
    gNonlData->SetLineColor(kBlue+1);
    gNonlData->SetLineWidth(3);
    gNonlData->SetMarkerSize(1.2);
    gNonlCalc->SetMarkerStyle(21);
    gNonlCalc->SetMarkerColor(kRed+1);
    gNonlCalc->SetMarkerSize(1.2);
    gNonlCalc->SetLineColor(kRed+1);
    gNonlCalc->SetLineWidth(3);

    TCanvas* c1 = new TCanvas("Nonlinearity", "Nonlinearity");
    c1->cd(); c1->SetGrid();
    gNonlData->SetTitle("Nonlinearity Fitting; Etrue/MeV; Nonlinearity");
    //gNonlData->GetYaxis()->SetRangeUser(0.01,0.045);
    gNonlData->Draw("APL");
    gNonlCalc->Draw("P SAME");
    TLegend* led = new TLegend();
    led->SetFillColor(kWhite);
    led->AddEntry(gNonlData, "data", "PL");
    led->AddEntry(gNonlCalc, "calc", "PL");
    led->Draw("SAME");

    c1->SaveAs("GamNLFit.root");    
}
