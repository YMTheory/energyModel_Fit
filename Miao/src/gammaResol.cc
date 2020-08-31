#include "gammaResol.hh"
#include "electronResol.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "junoParameters.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <TH1.h>
#include <TF1.h>
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

    m_loadData   = false;

    m_FitNL      = true;
    m_FitRes     = false;

    h_totalPE = new TH1D(m_name.c_str(), "", 10000, 0, 13000);
}

gammaResol::~gammaResol(){
    delete h_totalPE;
}


void gammaResol::LoadElecNLData()
{
    ifstream in; in.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/MC_NonL.txt");

    string line; double tmp_nonl, tmp_E; int num = 0;
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> tmp_E >>tmp_nonl;
        elecEtrue[num] = tmp_E;
        elecNonl[num] = tmp_nonl*1481.06/m_scale; 
        num++;
    }
    in.close();

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
}


void gammaResol::LoadGammaNLData()
{
    //cout << " >>> Loading Naked Gamma Data <<< " << endl;

    ifstream in;
    in.open(junoParameters::gammaLSNL_File);
    string line;

    int num = 0;
    string tmp_name; double tmp_E, tmp_totPE, tmp_totPESigma, tmp_EvisError;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_name >> tmp_E >> tmp_totPE >> tmp_totPESigma >> tmp_EvisError ;
        if(tmp_name == m_name) {
            m_Etrue = tmp_E;
            m_nonlData = tmp_totPE/m_scale/tmp_E;
            m_nonlDataErr = tmp_EvisError*tmp_totPE/m_scale/tmp_E;
            m_Evis = tmp_totPE/m_scale;
            m_resData = tmp_totPESigma/tmp_totPE;
            m_resDataErr = 0.01 * tmp_totPESigma/tmp_totPE;
        }
    }
    in.close();
    
}


void gammaResol::LoadData()
{
    LoadElecNLData    ();
    LoadGammaNLData   ();
    //LoadResData       ();

    m_loadData = true;
}




void gammaResol::calcGammaNPE()
{
    if(!m_loadData) LoadData();
    
    string filename = "./data/naked_gamma/"+m_name+"_all.txt";
    ifstream in; in.open(filename.c_str());
    if(!in) {cout << "No such file : " << filename << endl;}
    //cout << " >>> Reading " << filename << endl;
    Double_t tmp_elec; Int_t index = 0;
    string line; int num = 0;
    while(getline(in,line)) {
        num = 0;
        istringstream ss(line);
        while(ss>>tmp_elec) {
            EprmElec.push_back(tmp_elec); num++;
        }
        double tmp_Evis = 0;  double tmp_Esigma = 0.; double tmp_Etrue = 0;
        for(int j=0; j<num; j++) {
            int idx ;
            if(EprmElec[j]<0.01) { idx = int(EprmElec[j]/0.001)+1; }
            else {idx = int(EprmElec[j]/0.01)+10;}
            double elec_nonl;
            if(idx==0) {elec_nonl = elecNonl[idx];}
            else {elec_nonl = interpolate_nonl(idx, EprmElec[j]);}
            double elec_nonl1 = electronQuench::ScintillatorNL(EprmElec[j])+electronCerenkov::getCerenkovPE(EprmElec[j]);
            //cout <<"nonl: " << EprmElec[j] << " " << elec_nonl << " " << elec_nonl1 << endl;
            tmp_Etrue += EprmElec[j];
            tmp_Evis += elec_nonl1*EprmElec[j];
            tmp_Esigma += TMath::Power(electronResol::Resolution(EprmElec[j])*EprmElec[j],2);
            //tmp_Evis += EprmElec[j];
        }

        
        m_mean[index] = tmp_Evis*m_scale; 
        m_sigma[index] = TMath::Sqrt(tmp_Esigma)*m_scale;
        //cout << tmp_Etrue << " " << tmp_Evis << " " << m_mean[index] << endl;
        //cout << "prepare distributions: " << index << " " << m_mean[index] <<  " " << m_sigma[index] << endl;
        EprmElec.clear();
        index++; if(index==5000) break;

    }
    in.close();

    // 2-layer sample gamma totalPE distribution...
    double mean = 0; double sigma = 0;
    //double* m_sampleTotPE = new double[m_nSamples];
    for(int iSample=0; iSample<m_nSamples; iSample++) {
        int entry = int(gRandom->Uniform(0,5000));
        double totpe = gRandom->Gaus(m_mean[entry], m_sigma[entry]);
        //m_sampleTotPE[iSample] = totpe; 
        //mean += totpe;
        h_totalPE->Fill(totpe);
    } 
    h_totalPE->Fit("gaus", "Q");
    TF1* func = (TF1*)h_totalPE->GetFunction("gaus");
    //Esample = func->GetParameter(1);
    //tmp_sigma = func->GetParameter(2);
    m_resCalc = func->GetParameter(2)/func->GetParameter(1);
    m_nonlCalc = func->GetParameter(1)/m_scale/m_Etrue; 
    func->Delete();
    //cout << m_name << " " << m_Etrue << " "  << m_nonlCalc << " " << m_resCalc << " "  << endl;
    //mean /= m_nSamples;
    //for(int iSample=0; iSample<m_nSamples; iSample++) {
    //    sigma += TMath::Power((m_sampleTotPE[iSample]-mean), 2);
    //    h_totalPE->Fill(m_sampleTotPE[iSample]);
    //} sigma = TMath::Sqrt(sigma/(m_nSamples-1));

    //m_resCalc = sigma / mean;
    //m_nonlCalc = mean / m_scale/m_Etrue ;
    //delete []m_sampleTotPE;
}


double gammaResol::GetChi2()
{
    double chi2 = 0;
    
    // calculate totpe sigma
    calcGammaNPE();

    if(m_FitNL) {
        chi2 += TMath::Power( (m_nonlCalc - m_nonlData)/m_nonlDataErr, 2);
    }
    if(m_FitRes) {
        chi2 += TMath::Power( (m_resCalc - m_resData)/m_resDataErr, 2);
    }
    //cout << m_name << " chi2: " << chi2 << endl;
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


double gammaResol::interpolate_nonl(int idx, double E) {
    double delta = (elecNonl[idx]-elecNonl[idx-1])*(E-elecEtrue[idx-1]/1000)/(elecEtrue[idx]/1000-elecEtrue[idx-1]/1000) ;
    return delta+elecNonl[idx-1];
}




void gammaResol::check_nonl()
{
    LoadElecNLData();
    for(int i=0; i<800; i++) {
        double etrue = elecEtrue[i]/1000.;
        cout << etrue << " "<< elecNonl[i] << " " << electronQuench::ScintillatorNL(etrue)+electronCerenkov::getCerenkovPE(etrue) << endl; 
    } 
}
