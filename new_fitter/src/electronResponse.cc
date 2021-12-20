#include "electronResponse.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "junoParameters.hh"

#include <TFile.h>
#include <TGraph.h>

double electronResponse::m_SimEtrue[m_nSimData];
double electronResponse::m_SimNonl[m_nSimData];
double electronResponse::m_scale = 1503.664;
bool electronResponse::m_loadSimFile = false;
bool electronResponse::m_doFit = false;
bool electronResponse::m_loadResol = false;

TGraph* electronResponse::gSimData;
TF1* electronResponse::fElecNonl;
double electronResponse::m_p0 = 1.025;
double electronResponse::m_p1 = 0.1122;
double electronResponse::m_p2 = 1.394;
double electronResponse::m_p3 = 5.55e-4;

double electronResponse::m_q0 = 22.333;
double electronResponse::m_q1 = 0.982622;
double electronResponse::m_q2 = 6.12389e-05;

double electronResponse::m_ra = -2.17203e+00;
double electronResponse::m_rb = 1.31498e+03;
double electronResponse::m_rc = 1.60508e+02;

double electronResponse::m_c0 = 0 ; //-158.41;
double electronResponse::m_c1 = 1.77; //4.333;
double electronResponse::m_c2 = 5.77e-2; //0.00181;
double electronResponse::m_s0 = 1.;
double electronResponse::m_d0 = 0; //0;
double electronResponse::m_d1 = 0.20708522; //4.15105e-02;
double electronResponse::m_d2 = 0.00237574; //6.44493e-06;

double electronResponse::m_ma = 0.968;
double electronResponse::m_mb = 8.08e-3;
double electronResponse::m_mc = 0;

double electronResponse::m_na = 9.00687e-01;
double electronResponse::m_nb = 1.15790e-01;
double electronResponse::m_nc = 1.42189e+00;

double electronResponse::m_na1 = 6.13396e-01;
double electronResponse::m_nc1 = 1.15458e+00;

double electronResponse::m_n1 = 0.0035;
double electronResponse::m_n2 = 0.0057;

TGraphErrors* electronResponse::gMinElecNonl;
TGraphErrors* electronResponse::gMinElecNPE;
TGraphErrors* electronResponse::gElecResol;

double gElecResolFunc(double* x, double* p) {
    double E = x[0];
    double p0 = p[0];
    double p1 = p[1];
    double p2 = p[2];

    double sigma2 = p0 + p1*E + p2*E*E;
    if (sigma2<0)
        return 0;
    else 
        return TMath::Sqrt(sigma2);
}

double gNPESigmaFunc(double* x, double* p){
    double pe = x[0];
    double p0 = p[0];
    double p1 = p[1];
    double p2 = p[2];

    double sigma2 = p0 + p1*pe + p2*pe*pe;
    if (sigma2<0)
        return 0;
    else 
        return TMath::Sqrt(sigma2);
}


double gEvisSigma(double* x, double* p) {
    double pe = x[0];
    double ma = p[0];
    double mb = p[1];
    double mc = p[2];


    double sigma2 = mc*mc + ma*ma*pe + mb*mb*pe*pe;
    if (sigma2<0)
        return 0;
    else 
        return TMath::Sqrt(sigma2);
}


double gCerPESigma(double* x, double* p) {
    double E = x[0];
    double p0 = p[0]*p[0];
    double p1 = p[1]*p[1];
    double p2 = p[2]*p[2];

    double sigma2 = p0 + p1*E + p2*E*E;
    if (sigma2<0)
        return 0;
    else 
        return sigma2;
}

double gSctPESigma(double* x, double* p) {
    double E = x[0];
    double p0 = p[0];

    double sigma2 = p0*E ;
    if (sigma2<0)
        return 0;
    else 
        return sigma2;
}

double gNtotCov(double* x, double* p) {
    double N = x[0];
    double p0 = p[0]*p[0];
    double p1 = p[1]*p[1];
    double p2 = p[2]*p[2];

    double sigma2 = p0 + p1*N + p2*N*N;
    if (sigma2<0)
        return 0;
    else 
        return sigma2;
}

double gEvisCorr(double* x, double *p) {
    double E = x[0];
    double p0 = p[0];
    double p1 = p[1];
    return p0*E+p1;

}

double gEvisNew(double* x, double* p) {
    double E  = x[0];
    double p0 = p[0]*p[0];
    double p1 = p[1]*p[1];
    double n  = p[2];

    return TMath::Sqrt(p0*E + p1*TMath::Power(E, n));
}


double gEvisNew1(double* x, double* p) {
    double E  = x[0];
    double p0 = p[0]*p[0];
    double n  = p[1];

    return TMath::Sqrt( p0*TMath::Power(E, n));
}

TF1* electronResponse::fElecResol = new TF1("fElecNonl", gElecResolFunc, 0, 16, 3);
TF1* electronResponse::fCerPESigma = new TF1("fCerPESigma", gCerPESigma, 0, 1600, 3);
TF1* electronResponse::fSctPESigma = new TF1("fSctPESigma", gSctPESigma, 0, 1600, 1);
TF1* electronResponse::fNPESigma = new TF1("fNPESigma", gNPESigmaFunc, 0, 26000, 3);
TF1* electronResponse::fEvisSigma = new TF1("fEvisSigma", gEvisSigma, 0, 100000, 3);
TF1* electronResponse::fNtotCov = new TF1("fNtotCov", gNtotCov, 0, 26000, 3);
TF1* electronResponse::fEvisCorr = new TF1("fEvisCorr", gEvisCorr, 0, 16, 2);
TF1* electronResponse::fEvisNew  = new TF1("fEvisNew", gEvisNew, 0, 100000, 3);
TF1* electronResponse::fEvisNew1  = new TF1("fEvisNew1", gEvisNew1, 0, 100000, 2);

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

void electronResponse::loadMinSimElecNonl()
{
    gMinElecNonl = new TGraphErrors();
    gMinElecNPE  = new TGraphErrors();
    ifstream in;
    in.open("./data/electron/electron_response.txt");
    if (!in) std::cout << " >>> No electron response file!!! <<< " << std::endl;
    string line;
    double Etrue, totpe, totpe_err, sigma, sigma_err;
    int index = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> Etrue >> totpe >> totpe_err >> sigma >> sigma_err;
        gMinElecNonl->SetPoint(index, Etrue, totpe/Etrue/m_scale);
        gMinElecNonl->SetPointError(index, 0, 0.001);
        gMinElecNPE->SetPoint(index, Etrue, totpe);
        gMinElecNPE->SetPointError(index, 0, totpe_err);
        index++;
    }

    in.close();
}






void electronResponse::FuncConstruct()
{
    fElecNonl = new TF1("fElecNonl", "([0]+[3]*x)/(1+[1]*TMath::Exp(-[2]*x))", 0.01, 8);
    SetParameters();
    cout << "Empirical Electron Response Function has been constructed..." << endl;
}

void electronResponse::SetParameters()
{
    //fElecNonl->SetParameter(0, m_p0);
    //fElecNonl->SetParameter(1, m_p1);
    //fElecNonl->SetParameter(2, m_p2);
    //fElecNonl->SetParameter(3, m_p3);

    fElecResol->SetParameter(0, m_ra);
    fElecResol->SetParameter(1, m_rb);
    fElecResol->SetParameter(2, m_rc);

    fCerPESigma->SetParameter(0, m_c0);
    fCerPESigma->SetParameter(1, m_c1);
    fCerPESigma->SetParameter(2, m_c2);

    fSctPESigma->SetParameter(0, m_s0);


    fNPESigma->SetParameter(0, m_q0);
    fNPESigma->SetParameter(1, m_q1);
    fNPESigma->SetParameter(2, m_q2);

    fEvisSigma->SetParameter(0, m_ma);
    fEvisSigma->SetParameter(1, m_mb);
    fEvisSigma->SetParameter(2, m_mc);

    fNtotCov->SetParameter(0, m_d0); 
    fNtotCov->SetParameter(1, m_d1); 
    fNtotCov->SetParameter(2, m_d2); 

    fEvisNew->SetParameter(0, m_na);
    fEvisNew->SetParameter(1, m_nb);
    fEvisNew->SetParameter(2, m_nc);


    fEvisNew1->SetParameter(0, m_na1);
    fEvisNew1->SetParameter(1, m_nc1);


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


void electronResponse::FitPlot()
{
    TGraph* gNom = new TGraph();
    TGraph* gFit = new TGraph();
    gNom->SetName("nom");
    gFit->SetName("fit");


    electronQuench::setkA(1);
    electronQuench::setEnergyScale(3300.371/2.223);
    electronCerenkov::setkC(1);
    electronCerenkov::setEnergyScale(3300.371/2.223);
    for(int i=0; i<500; i++) {
        double Etrue = 10./500 * (i+1);
        double nonl_nom = electronResponse::getElecNonl(Etrue);
        gNom->SetPoint(i, Etrue, nonl_nom);
    }

    gNom->SetLineWidth(2);
    gNom->SetLineColor(kBlue+1);


    electronQuench::setkA(1.01959e+00);
    electronQuench::setEnergyScale(3300.371/2.223);
    electronCerenkov::setkC(7.09311e-01);
    electronCerenkov::setEnergyScale(3300.371/2.223);
    for(int i=0; i<500; i++) {
        double Etrue = 10./500 * (i+1);
        double nonl_fit = electronResponse::getElecNonl(Etrue);
        gFit->SetPoint(i, Etrue, nonl_fit);
    }

    gFit->SetLineWidth(2);
    gFit->SetLineColor(kGreen+1);


    TFile* out = new TFile("ElecNonlCompare.root", "recreate");
    gNom->Write();
    gFit->Write();
    out->Close();


}

void electronResponse::loadElecResol()
{
    gElecResol = new TGraphErrors();
    ifstream in; in.open(junoParameters::electronResol_File.c_str());
    string line;
    double tmp_E, tmp_mu, tmp_muerr, tmp_sigma, tmp_sigmaerr, tmp_resol, tmp_resolerr;
    int index = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmp_E >> tmp_mu >> tmp_muerr >> tmp_sigma >> tmp_sigmaerr >> tmp_resol >> tmp_resolerr;
        gElecResol->SetPoint(index, tmp_mu, tmp_sigma);
        index++;
    }
    in.close();

    m_loadResol = true;
}





void electronResponse::SaveElecNonl(double* par)
{
    int iPar = 0;
    electronQuench::setEnergyScale(par[iPar]);                  iPar++;
    electronQuench::setBirk1(par[iPar]);                        iPar++;
    electronCerenkov::setkC(par[iPar]);                         iPar++;
    
    TFile* ff = new TFile("elec_nonl_onlyGam.root", "recreate");
    TGraph* gg = new TGraph();
    gg->SetName("nonl");
    for (int i=0; i<1490; i++) {
        double E = 0.01*(i+1);
        double npe = electronQuench::ScintillatorPE(E) + electronCerenkov::getCerPE(E);
        gg->SetPoint(i, E, npe);
    }
    gg->Write();
    ff->Close();
}




