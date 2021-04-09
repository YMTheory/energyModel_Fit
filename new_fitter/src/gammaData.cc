#include "gammaData.hh"
#include "junoParameters.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>

using namespace std;

//std::string gammaData::m_calcOption = "prmelec";
std::string gammaData::m_calcOption = junoParameters::m_calcOption;

//std::string gammaData::m_nonlMode = "histogram";
std::string gammaData::m_nonlMode = junoParameters::m_nonlMode;

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

    m_calcOption = junoParameters::m_calcOption;
    m_nonlMode   = junoParameters::m_nonlMode;

    m_loadData = false;

    hEmean = new TH1D(m_name.c_str(), "", 10000, 0, 10);
    hPEmean = new TH1D((m_name+"pe").c_str(), "", 10000, 0, 20000);

}

gammaData::~gammaData() 
{}

void gammaData::LoadGammaData()
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
            m_nonlDataErr = 0.001;
            m_Evis = tmp_totPE/scale;
            m_resData = tmp_totPESigma/tmp_totPE;
            //m_resDataErr = 0.01 * tmp_totPESigma/tmp_totPE;
            m_resDataErr = TMath::Sqrt(tmp_totPESigmaerr*tmp_totPESigmaerr/tmp_totPE/tmp_totPE + tmp_totPEerr*tmp_totPEerr*tmp_totPESigma*tmp_totPESigma/tmp_totPE/tmp_totPE/tmp_totPE/tmp_totPE);
        }
    }
    in.close();
}

void gammaData::LoadPrimElecDist()
{
    cout << " >>> Load Primary Electron Distribution <<< " << endl;

    TFile* file = new TFile(junoParameters::gammaPdf_File.c_str(), "read");
    if (!file) cout << " No such input primary electron file " << endl;
    //string pdfName = m_name;
    string pdfName = "gamma" + m_name;
    TH1D* gGammaPdf = (TH1D*)file->Get(pdfName.c_str());
    //double intvalue = gGammaPdf->Integral() ;
    //gGammaPdf->Scale(1/intvalue);
    if (!gGammaPdf) cout << " No such Pdf : " << pdfName << endl;

    for(int i=0; i<gGammaPdf->GetNbinsX(); i++) {
        m_pdf_eTrue[i] = gGammaPdf->GetBinCenter(i+1);
        m_pdf_prob[i]  = gGammaPdf->GetBinContent(i+1);
        if (m_pdf_prob[i] == 0 ) { 
            m_max_eTrue = i;
            break;
        }
    }

    file->Close();

}

void gammaData::LoadPrimElecSamples()
{
    cout << " >>> Load Primary Electron in Single Event <<< " << endl;
    string filename = "./data/gamma/" + m_name + "_all.root";
    TFile* file = new TFile(filename.c_str(), "read");
    if (!file) cout << " No such input file: " << filename << endl;
    elec_hist = (TH2D*)file->Get(m_name.c_str());

}


void gammaData::LoadData()
{
    LoadGammaData();

    if (m_calcOption == "prmelec")
        LoadPrimElecDist();

    if (m_calcOption == "twolayer")
        LoadPrimElecSamples();

    m_loadData = true;
}


void gammaData::calcGammaResponse()
{
    if (!m_loadData) LoadData();
    m_scale = electronQuench::getEnergyScale();

    if (m_calcOption == "prmelec") {
        double numerator = 0;
        double denominator = 0;
        
        double numerator1 = 0;
        double denominator1 = 0;
        for(int iBin=1; iBin<m_max_eTrue; iBin++) {
            double E1 = m_pdf_eTrue[iBin-1];
            double E2 = m_pdf_eTrue[iBin];
            double prob1 = m_pdf_prob[iBin-1];
            double prob2 = m_pdf_prob[iBin];

            // based on DYB method...
            double fNL1, fNL2;

            if(m_nonlMode == "histogram") {
                fNL1 = electronQuench::ScintillatorNL(E1) + electronCerenkov::getCerenkovPE(E1);
                fNL2 = electronQuench::ScintillatorNL(E2) + electronCerenkov::getCerenkovPE(E2);
            }

            if (m_nonlMode == "analytic") {
                fNL1 = electronResponse::calcElecNonl(E1);
                fNL2 = electronResponse::calcElecNonl(E2);
            }

            numerator += (prob1*E1*fNL1 + prob2*E2*fNL2) * (E2-E1) / 2.;
            denominator += (prob1*E1 + prob2*E2) * (E2-E1) / 2.;

            double NPE1, NPE2;
            if (m_nonlMode == "histogram" ) {
                NPE1 = electronQuench::ScintillatorPE(E1) + electronCerenkov::getCerPE(E1);
                NPE2 = electronQuench::ScintillatorPE(E2) + electronCerenkov::getCerPE(E2);
            }

            cout << "mid-process " << (prob1*NPE1 + prob2*NPE2) *(E2-E1)/2 << " " << (prob1 + prob2) * (E2-E1)/2 << endl;
            numerator1 += (prob1*NPE1 + prob2*NPE2) * (E2-E1) /2;
            denominator1 += (prob1 + prob2) *(E2-E1) /2;
        }

        if (denominator == 0) cout << "Errors Happend While Using GammaPdf Calculation..." << endl;

        m_nonlCalc = numerator / denominator;
        m_totpeCalc = m_nonlCalc * m_Etrue * m_scale;
        m_nonlCalc1 = numerator1 /denominator1 / m_scale / m_Etrue;
        cout << numerator1 << " " << denominator1 << " " << m_scale << " " << m_Etrue << endl;

    }  else if (m_calcOption == "twolayer") {
        for (int iSample=0; iSample<m_nSamples; iSample++) {
            hEmean->Reset();
            hPEmean->Reset();
        
            // apply Nonlinearity curve
            double tmp_mean = 0;
            double tmp_pe= 0;
            for (int iSec=0; iSec<100; iSec++) {
                double tmp_E = elec_hist->GetBinContent(iSample+1, iSec+1);
                if (tmp_E == 0) break;
                double tmp_Edep;
                double tmp_sigpe;
                if(m_nonlMode == "histogram") {
                    tmp_Edep = tmp_E * electronResponse::getElecNonl(tmp_E);
                    tmp_sigpe = (electronQuench::ScintillatorPE(tmp_E) +electronCerenkov::getCerPE(tmp_E)); // alreadly include nonl
                }
                if(m_nonlMode == "analytic")
                    tmp_Edep = tmp_E * electronResponse::calcElecNonl(tmp_E);
                tmp_mean += tmp_Edep;
                tmp_pe += tmp_sigpe;
                
            }
            m_mean[iSample] = tmp_mean;
            //m_meanpe[iSample] = tmp_pe;
            
            hEmean->Fill(tmp_mean);
            hPEmean->Fill(tmp_pe);
        }

        // calculate pe distribution
        double mean_Edep = 0;
        double mean_pe = 0;
        for (int i=0; i<m_nSamples; i++) {
            mean_Edep += m_mean[i];
            //mean_pe += m_meanpe[i];
        }
        mean_Edep /= m_nSamples;
        //mean_pe /= m_nSamples;
        //m_totpeCalc = mean_pe;
        m_nonlCalc = mean_Edep / m_Etrue;
        m_nonlCalc1 = hEmean->GetMean() / m_Etrue;
        m_totpeCalc = hPEmean->GetMean();
        //m_nonlCalc = hEmean->GetMean() / m_Etrue;
        //m_nonlCalc1 = m_totpeCalc / m_scale / m_Etrue ;
    }

    // pull term for gamma :
    m_nuGamma = junoParameters::m_nuGamma;
    m_nonlCalc *= (1+m_nuGamma);
 
}


double gammaData::GetChi2()
{
    double chi2 = 0;

    // calculate totpe sigma
    calcGammaResponse();

    chi2 += (m_nonlCalc - m_nonlData) * (m_nonlCalc - m_nonlData) / m_nonlDataErr / m_nonlDataErr;

    return chi2;
}

