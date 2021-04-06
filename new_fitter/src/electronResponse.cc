#include "electronResponse.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

#include <TFile.h>

double electronResponse::m_SimEtrue[m_nSimData];
double electronResponse::m_SimNonl[m_nSimData];
double electronResponse::m_scale = 3300.371/2.223;
bool electronResponse::m_loadSimFile = false;
bool electronResponse::m_doFit = false;

TGraph* electronResponse::gSimData;
TF1* electronResponse::fElecNonl;
double electronResponse::m_p0 = 1.025;
double electronResponse::m_p1 = 0.1122;
double electronResponse::m_p2 = 1.394;
double electronResponse::m_p3 = 5.55e-4;

double electronResponse::getElecNonl(double Etrue)
{
    return electronQuench::ScintillatorNL(Etrue) + electronCerenkov::getCerenkovPE(Etrue);
}


void electronResponse::loadSimElecNonl()
{
    gSimData = new TGraph();
    ifstream in; 
    in.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/totPE_smearing.txt");
    string line;
    double Etrue, nonl, totpe, totpe_sigma;
    int index = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> Etrue >> totpe >> totpe_sigma;
        m_SimEtrue[index] = Etrue/1000.;
        if (Etrue == 0)
            m_SimNonl[index] = 0;
        else{
            nonl = totpe/Etrue/m_scale*1000.;
            m_SimNonl[index] = nonl;
        }
        gSimData->SetPoint(index, m_SimEtrue[index], m_SimNonl[index]);
        index++;
    }
    in.close();

    m_loadSimFile = true;
    return;
}

void electronResponse::FuncConstruct()
{
    fElecNonl = new TF1("fElecNonl", "([0]+[3]*x)/(1+[1]*TMath::Exp(-[2]*x))", 0.01, 8);
    SetParameters();
    cout << "Empirical ElecNonl Function has been constructed..." << endl;
}

void electronResponse::SetParameters()
{
    fElecNonl->SetParameter(0, m_p0);
    fElecNonl->SetParameter(1, m_p1);
    fElecNonl->SetParameter(2, m_p2);
    fElecNonl->SetParameter(3, m_p3);
}


void electronResponse::Plot()
{
    cout << ">>> Draw Analytic Electron Nonlinearity Curve <<< "<< endl;
    if (not m_loadSimFile) loadSimElecNonl();
    //if (not m_doFit) EmpiricalFit();

    gSimData->SetName("elec");
    gSimData->SetLineColor(kBlue+1);
    gSimData->SetLineWidth(2);

    fElecNonl->SetLineColor(kOrange+1);
    fElecNonl->SetLineWidth(2);

    TFile* out = new TFile("simElecNonl.root", "recreate");
    gSimData->Write();
    fElecNonl->Write();
    out->Close();
}


void electronResponse::EmpiricalFit()
{
    if(not m_loadSimFile) loadSimElecNonl();
    fElecNonl->SetParameters(1.025, 0.1122, 1.394, 5.55e-4);
    gSimData->Fit(fElecNonl, "R");

    m_doFit = true;
}
