#include "gammaData.hh"
#include "junoParameters.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include <TFile.h>
#include <TH1.h>

using namespace std;

double gammaData::m_pdf_eTrue[m_nMaxPdf];
double gammaData::m_pdf_prob[m_nMaxPdf];
double gammaData::m_max_eTrue;

string gammaData::m_calcOption = junoParameters::gammaNLOption;

gammaData::gammaData( std::string name,
                        double minPE,
                        double maxPE,
                        int nbins
                      ) 
{
    m_name = name;
    m_minPE = minPE;
    m_maxPE = maxPE;
    m_nbins = nbins;

    m_loadData = false;
}

gammaData::~gammaData() {

}

void gammaData::LoadGammaData()
{
    cout << " >>> Loading Naked Gamma " << m_name << " Data" << endl;

    ifstream in; in.open(junoParameters::gammaLSNL_File);
    string line;

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

void gammaData::LoadPrimElecDist()
{
    cout << " >>> Load Primary Electron Distribution <<< " << endl;

    TFile* file = new TFile(junoParameters::gammaPdf_File.c_str(), "read");
    if (!file) cout << " No such input primary electron file " << endl;
    string pdfName = "gamma" + m_name;
    TH1D* gGammaPdf = (TH1D*)file->Get(pdfName.c_str());
    if (!gGammaPdf) cout << " No such Pdf : " << pdfName << endl;

    for(int i=0; i<gGammaPdf->GetNbinsX(); i++) {
        m_pdf_eTrue[i] = gGammaPdf->GetBinCenter(i+1);
        m_pdf_prob[i]  = gGammaPdf->GetBinContent(i+1);
        if (m_pdf_prob[i] > 0) m_max_eTrue = i;
    }

    file->Close();
}


void gammaData::LoadData()
{
    LoadGammaData();

    if (m_calcOption == "prmelec")
        LoadPrimElecDist();

    m_loadData = true;
}


void gammaData::calcGammaResponse()
{
    if (!m_loadData) LoadData();

    if (m_calcOption == "primelec") {
        double numerator = 0;
        double denominator = 0;
        
        for(int iBin=1; iBin<m_max_eTrue; iBin++) {
            double E1 = m_pdf_eTrue[iBin-1];
            double E2 = m_pdf_eTrue[iBin];
            double prob1 = m_pdf_prob[iBin-1];
            double prob2 = m_pdf_prob[iBin];

            double fNL1 = electronQuench::ScintillatorNL(E1) + electronCerenkov::getCerenkovPE(E1);
            double fNL2 = electronQuench::ScintillatorNL(E2) + electronCerenkov::getCerenkovPE(E2);

            numerator += (prob1*E1*fNL1 + prob2*E2*fNL2) * (E2-E1) / 2.;
            denominator += (prob1*E1 + prob2*E2) * (E2-E1) / 2.;
        }

        if (denominator == 0) cout << "Errors Happend While Using GammaPdf Calculation..." << endl;

        m_nonlCalc = numerator / denominator;
    }
}


double gammaData::GetChi2()
{
    double chi2 = 0;

    // calculate totpe sigma
    calcGammaResponse();

    if (m_calcOption == "primelec") {
        chi2 += (m_nonlCalc - m_nonlData) * (m_nonlCalc - m_nonlData) / m_nonlDataErr / m_nonlDataErr;
    }
    
    return chi2;
}

