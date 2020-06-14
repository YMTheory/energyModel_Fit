#include "gammaResol.hh"
#include "electronResol.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <TH1.h>
#include <TFile.h>
#include <TRandom.h>
#include <TMath.h>

using namespace std;

gammaResol::gammaResol(std::string name,
                       double minPE,
                       double maxPE,
                        int nbins) {
    m_name = name;
    m_minPE = minPE;
    m_maxPE = maxPE;
    m_nbins = nbins;

    m_loadNL    = false;
    m_loadRes   = false;

    h_totalPE = new TH1D(m_name.c_str(), "", m_nbins, m_minPE, m_maxPE);
}

gammaResol::~gammaResol(){
    delete h_totalPE;
}


void gammaResol::LoadElecNLData()
{
    ifstream in; in.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/MC_NonL.txt");
    string line; double tmp_nonl; int num = 0;
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> tmp_nonl;
        elecNonl[num] = tmp_nonl; num++;
    }
    in.close();

    m_loadNL = true;
}


void gammaResol::LoadResData()
{
    ifstream in; in.open("./data/naked_gamma/gamma_resol.txt");
    if(!in) {cout << "No Such A File" << endl; }
    string line;
    string tmp_name; double tmp_E; double tmp_res; double tmp_resErr;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmp_name >> tmp_E >> tmp_res >> tmp_resErr;
        if (tmp_name == m_name) {
            m_Evis = tmp_E; m_resData = tmp_res; m_resDataErr = tmp_resErr; break;
        }
    }

    in.close();
    m_loadRes = true;
}


void gammaResol::calcGammaNPE()
{
    if(!m_loadNL) LoadElecNLData();
    
    string filename = "./data/naked_gamma/"+m_name+".txt";
    ifstream in; in.open(filename.c_str());
    if(!in) {cout << "No such file : " << filename << endl;}
    Double_t tmp_elec; Int_t index = 0;
    string line; int num = 0;
    while(getline(in,line)) {
        num = 0;
        istringstream ss(line);
        while(ss>>tmp_elec) {
            EprmElec.push_back(tmp_elec); num++;
        }
        double tmp_Evis = 0;  double tmp_Esigma = 0.;
        for(int j=0; j<num; j++) {
            int idx = int(EprmElec[j]/m_NLResol);
            double elec_nonl = elecNonl[idx];
            tmp_Evis += elec_nonl*EprmElec[j];
            tmp_Esigma += TMath::Power(electronResol::Resolution(EprmElec[j])*EprmElec[j],2);
            //tmp_Evis += EprmElec[j];
        }

        //Evis_calc[index] = tmp_Evis; index++;
        m_mean[index] = tmp_Evis*m_scale; 
        m_sigma[index] = TMath::Sqrt(tmp_Esigma)*m_scale;
        index++;
        EprmElec.clear();

    }
    in.close();

    double mean = 0; double sigma = 0;
    double* m_sampleTotPE = new double[m_nSamples];
    for(int iSample=0; iSample<m_nSamples; iSample++) {
        int entry = int(gRandom->Uniform(0,1000));
        double totpe = gRandom->Gaus(m_mean[entry], m_sigma[entry]);
        m_sampleTotPE[iSample] = totpe;
        mean += totpe;
    } mean /= m_nSamples;
    for(int iSample=0; iSample<m_nSamples; iSample++) {
        sigma += TMath::Power((m_sampleTotPE[iSample]-mean), 2);
        h_totalPE->Fill(m_sampleTotPE[iSample]);
    } sigma = TMath::Sqrt(sigma/(m_nSamples-1));

    m_resCalc = sigma / mean;
    delete []m_sampleTotPE;
}


//void gammaResol::predGammaNPE()
//{
//    string filename = "./data/naked_gamma/singleEvent_"+m_name+".txt";
//    ifstream in;
//    in.open(filename.c_str());
//    if(!in) {cout << "No Such singleEvent File: "<< filename << endl;}
//    double tmp_mean; double tmp_sigma;
//    string line; int index = 0;
//    while(getline(in, line)) {
//        istringstream ss(line);
//        ss >> tmp_mean >> tmp_sigma;
//        m_mean[index] = tmp_mean;
//        m_sigma[index] = tmp_sigma; index++;
//    }
//    in.close();
//    
//    
//    double mean = 0; double sigma = 0;
//    for(int iSample=0; iSample<m_nSamples; iSample++) {
//        int entry = int(gRandom->Uniform(0,1000));
//        double totpe = gRandom->Gaus(m_mean[entry], m_sigma[entry]);
//        m_sampleTotPE[iSample] = totpe;
//        mean += totpe;
//    } mean /= m_nSamples;
//    for(int iSample=0; iSample<m_nSamples; iSample++) {
//        sigma += TMath::Power((m_sampleTotPE[iSample]-mean), 2);
//        h_totalPE->Fill(m_sampleTotPE[idx]);
//    } sigma = TMath::Sqrt(sigma/(m_nSamples-1));
//
//    m_resCalc = sigma / mean;
//    
//}

double gammaResol::GetChi2()
{
    double chi2 = 0;
    
    // calculate totpe sigma
    calcGammaNPE();
    if(!m_loadRes) LoadResData();

    chi2 = TMath::Power( (m_resCalc - m_resData)/m_resDataErr, 2);
    //cout << "mean totpe: " << mean << ", sigma: " << sigma << endl;
    return chi2;
}


void gammaResol::Plot()
{

    //h_totalPE = new TH1D(m_name.c_str(), "", m_nbins, m_minPE, m_maxPE);
    //for(int idx=0; idx<m_nSamples; idx++) {
    //    h_totalPE->Fill(m_sampleTotPE[idx]);
    //}

    h_totalPE->SetLineColor(kBlue+1);
    h_totalPE->SetMarkerColor(kBlue+1);
    h_totalPE->SetMarkerStyle(20);
    h_totalPE->SetMarkerSize(0.6);
    //h_totalPE->Draw("PEX0");

    TFile* file = new TFile("./output/gamma/predCs137.root", "recreate");
    h_totalPE->Write();
    file->Close();
}
