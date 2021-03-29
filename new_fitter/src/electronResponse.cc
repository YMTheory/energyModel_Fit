#include "electronResponse.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

#include <TFile.h>

double electronResponse::m_SimEtrue[m_nSimData];
double electronResponse::m_SimNonl[m_nSimData];
double electronResponse::m_scale = 3350/2.22;
bool electronResponse::m_loadSimFile = false;
bool electronResponse::m_doFit = false;

TGraph* electronResponse::gSimData;
TF1* electronResponse::func;

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

void electronResponse::Plot()
{
    if (not m_loadSimFile) loadSimElecNonl();
    if (not m_doFit) EmpiricalFit();

    gSimData->SetName("elec");
    gSimData->SetLineColor(kBlue+1);
    gSimData->SetLineWidth(2);

    func->SetLineColor(kOrange+1);
    func->SetLineWidth(2);

    TFile* out = new TFile("simElecNonl.root", "recreate");
    gSimData->Write();
    func->Write();
    out->Close();
}


void electronResponse::EmpiricalFit()
{
    if(not m_loadSimFile) loadSimElecNonl();
    func = new TF1("lsnl_func", "([0]+[3]*x)/(1+[1]*TMath::Exp(-[2]*x))", 0.01, 8);
    func->SetParameters(1.025, 0.1122, 1.394, 5.55e-4);
    gSimData->Fit(func, "R");

    m_doFit = true;
}