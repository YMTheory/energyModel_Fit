#include "gammaResponse.hh"
#include "junoParameters.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TFile.h>
#include <TRandom.h>
#include <TF1.h>
#include <TTree.h>

gammaResponse::gammaResponse(string name, int nBins, double peMin, double peMax) {
    m_name = name;
    m_nBins = nBins;
    m_peMin = peMin;
    m_peMax = peMax;

    hCalc = new TH1D((m_name+"_calc").c_str(), "", m_nBins, m_peMin, m_peMax); 
    hData = new TH1D((m_name+"_data").c_str(), "", m_nBins, m_peMin, m_peMax); 

    m_loadData = false;
    m_loadPrm = false;
    m_doSpecFit = true;
}

gammaResponse::~gammaResponse()
{
    delete hPrmElec;
    delete hPrmPosi;

    delete hCalc;
}


void gammaResponse::LoadData()
{
    cout << " >>> Loading Naked Gamma " << m_name << " Data <<< " << endl;

    ifstream in; in.open(junoParameters::gammaLSNL_File);
    string line;

    double scale = 3300.371/2.223;
    string tmp_name; double tmp_E, tmp_totPE, tmp_totPESigma, tmp_EvisError, tmp_totPEerr, tmp_totPESigmaerr;
    while(getline(in,line)){
        istringstream ss(line);
        //ss >> tmp_name >> tmp_E >> tmp_totPE >> tmp_totPESigma >> tmp_EvisError ;
        ss >> tmp_name >> tmp_E >> tmp_totPE >> tmp_totPEerr >> tmp_totPESigma >> tmp_totPESigmaerr;
        if(tmp_name == m_name) {
            m_Etrue = tmp_E;
            m_totpeData = tmp_totPE;
            m_nonlData = tmp_totPE/scale/tmp_E;
            //m_nonlDataErr = tmp_EvisError*tmp_totPE/scale/tmp_E;
            //m_nonlDataErr = tmp_totPEerr/scale/tmp_E;
            m_nonlDataErr = 0.002;
            m_Evis = tmp_totPE/scale;
            m_resData = tmp_totPESigma/tmp_totPE;
            //m_resDataErr = 0.01 * tmp_totPESigma/tmp_totPE;
            m_resDataErr = TMath::Sqrt(tmp_totPESigmaerr*tmp_totPESigmaerr/tmp_totPE/tmp_totPE + tmp_totPEerr*tmp_totPEerr*tmp_totPESigma*tmp_totPESigma/tmp_totPE/tmp_totPE/tmp_totPE/tmp_totPE);
            break;
        }
    }
    in.close();

    if (m_doSpecFit) {
        TFile* infile = new TFile(("./data/gamma/spectrum/" + m_name + "_totpe.root").c_str(), "read");
        if(!infile) cout << "No such gamma spectrum file!" << endl;
        double m_totpe;
        TTree* tt = (TTree*)infile->Get(m_name.c_str());
        tt->SetBranchAddress("totpe", &m_totpe);
        for(int i=0; i<tt->GetEntries(); i++) {
            tt->GetEntry(i);
            hData->Fill(m_totpe);
        }
    }


    LoadPrmBeta();
}

void gammaResponse::LoadPrmBeta()
{
    cout << " >>> Load Primary Electron in Single Event <<< " << endl;
    string filename = "./data/gamma/" + m_name + "_new.root";
    TFile* file = new TFile(filename.c_str(), "read");
    if (!file) cout << " No such input file: " << filename << endl;
    hPrmElec = (TH2D*)file->Get((m_name+"_elec").c_str());
    hPrmPosi = (TH2D*)file->Get((m_name+"_posi").c_str());
    
    m_loadPrm = true;
}


double gammaResponse::SampleGamEnergy(int index)
{
    if (index >= m_nEvents or index < 0) { cout << "Incorrect Index !" << endl; return 0; } 
    else {
        double tmp_pe = 0;
        double tmp_sigma = 0;
        // electron
        for (int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmElec->GetBinContent(index+1, iSec+1);    
            if (tmp_E == 0) break;
            tmp_pe += electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E);
            tmp_sigma += TMath::Power(electronResponse::gElecResol->Eval(tmp_E), 2);
        }

        // positron
        for(int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmPosi->GetBinContent(index+1, iSec+1);
            if (tmp_E == 0) break;
            tmp_pe += electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E) + 695.53*2;
            tmp_sigma += TMath::Power(electronResponse::gElecResol->Eval(tmp_E), 2);
            tmp_sigma += 2*TMath::Power(27.513, 2);
            
        }

        tmp_sigma = TMath::Sqrt(tmp_sigma);
    
        double sample_pe = gRandom->Gaus(tmp_pe, tmp_sigma);
        //if (sample_pe < 8000) { // for small p.e. events check...
        //    cout << index << " " << tmp_pe << " " << tmp_sigma <<   endl;
        //}
        return sample_pe;
    }
}




void gammaResponse::calcGamResponse()
{
    if (not electronResponse::getLoadResol()) electronResponse::loadElecResol();
    hCalc->Reset();

    for (int iSample = 0; iSample<m_nSamples; iSample++) {
        int index = gRandom->Integer(5000);
        double sample_pe = SampleGamEnergy(index);
        hCalc->Fill(sample_pe); 
    }
    //hCalc->Fit("gaus", "Q0");
    //double pe_mean  = hCalc->GetFunction("gaus")->GetParameter(1);
    //double pe_sigma = hCalc->GetFunction("gaus")->GetParameter(2);

    double pe_mean  = hCalc->GetMean();
    double pe_sigma = hCalc->GetStdDev();

    m_nonlCalc = pe_mean / electronQuench::getEnergyScale() / m_Etrue;
    m_resCalc = pe_sigma / pe_mean;

}


double gammaResponse::GetChi2()
{
    double chi2 = 0;

    calcGamResponse();

    if (not m_doSpecFit) 
        chi2 += TMath::Power((m_nonlData - m_nonlCalc)/m_nonlDataErr, 2);

    else{
        for (int i=0; i<m_nBins; i++) {
            double m_err = hData->GetBinError(i);
            if (m_err != 0) {
                double m_data = hData->GetBinContent(i);
                double m_calc = hCalc->GetBinContent(i);

                chi2 += TMath::Power((m_calc - m_data)/m_err, 2);
                //cout << hData->GetBinCenter(i) << " " << hCalc->GetBinCenter(i) << " " << m_calc << " " << m_data << " " << m_err << " " << chi2 << endl;
            }
        }
    }
        
    return chi2 ;
}


void gammaResponse::SaveHist()
{
    string fileName = m_name + "hist.root";
    TFile* outfile = new TFile(fileName.c_str(), "recreate");
    hData->Write();
    hCalc->Write();
    outfile->Close();
}

