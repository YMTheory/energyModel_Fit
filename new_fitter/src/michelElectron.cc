#include "michelElectron.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

using namespace std;

michelElectron::michelElectron()
{

    m_dataEntries = 0;
    m_calcEntries = 0;

    m_nBin = 1000; //200;
    m_EdepMin = 0; //45;      // MeV
    m_EdepMax = 54;      // MeV
    m_binWidth = (m_EdepMax - m_EdepMin) / m_nBin;

    m_nBinData = 40;
    m_EvisMin = 50;      // MeV
    m_EvisMax = 57;      // MeV
    m_binDataWidth = (m_EvisMax - m_EvisMin) / m_nBinData;
    m_binCalcWidth = (m_EvisMax - m_EvisMin) / m_nBin;

    hEvisData = new TH1D("hEvisData", "hEvisData", m_nBinData, m_EvisMin, m_EvisMax);
    hEvisCalc = new TH1D("hEvisCalc", "hEvisCalc", m_nBinData, m_EvisMin, m_EvisMax);

    m_energyScale = 3134.078 / 2.223;

    hEdep = new TH1D("hEdep", "hEdep", m_nBin, m_EdepMin, m_EdepMax);
    hTmpEvisCalc = new TH1D("hTmpEvisCalc", "hTmpEvisCalc", m_nBin, m_EvisMin, m_EvisMax);

    gaus = new TF1("gaus", "1/(TMath::Sqrt(2*TMath::Pi())*[1]) * TMath::Exp(-(x-[0])*(x-[0])/2/[1]/[1])", -100, 100);
}

michelElectron::~michelElectron()
{
    delete hEdep;
    delete hEvisData;
    delete hEvisCalc;
    delete hTmpEvisCalc;

    delete gaus;
}


void michelElectron::Initialize()
{
    LoadMCTruth();
    LoadMCData();
}

void michelElectron::LoadMCTruth()
{
    TFile* ff = new TFile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/michel/michel_edep_LS.root", "read");
    //TFile* ff = new TFile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/michel/michel_totpe_LS_v5.root", "read");
    TTree* tt = (TTree*)ff->Get("michel");
    double m_edep; 
    tt->SetBranchAddress("edep", &m_edep);
    for(int i=0; i<tt->GetEntries(); i++) {
        tt->GetEntry(i);
        hEdep->Fill(m_edep);
    }

    cout << endl;
    cout << "********************************" << endl;
    cout << "        Load Michel MC Data       " << endl;
    cout << "********************************" << endl;
    cout << endl;

}


void michelElectron::LoadMCData()
{
    TFile* ff = new TFile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/michel/michel_totpe_LS_v7.root", "read");
    //TFile* ff = new TFile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/michel/michel_totpe_LS_v3.root", "read");
    TTree* tt = (TTree*)ff->Get("michel");
    int m_totpe;
    tt->SetBranchAddress("totpe", &m_totpe);
    for(int i=0; i<tt->GetEntries(); i++) {
        tt->GetEntry(i);
        hEvisData->Fill(m_totpe/m_energyScale);
        if(m_totpe/m_energyScale > m_EvisMin and m_totpe/m_energyScale < m_EvisMax) { m_dataEntries++;}
    }

}


void michelElectron::CalcVisibleEnergy()
{

    //electronQuench::setEnergyScale(1412.15);
    //electronQuench::setBirk1(6.431e-03);
    //electronCerenkov::setp0(0.906);
    //electronCerenkov::setp1(10.2);
    //electronCerenkov::setp2(0.017);
    //electronCerenkov::setp3(74.8);
    //electronCerenkov::setp4(-9.90);
    //electronResponse::setna(9.404e-01);
    //electronResponse::setnb(9.91e-02);
    //electronResponse::setnc(1.449);

    ////electronResponse::setma(1.019);
    ////electronResponse::setmb(3.3e-3);
    ////electronResponse::setmc(0);
    //electronResponse::SetParameters();

    // Reset:
    for(int i=0; i<m_nBin; i++) {
        hTmpEvisCalc->SetBinContent(i+1, 0);
    }
    for(int i=0; i<m_nBinData; i++) {
        hEvisCalc->SetBinContent(i+1, 0);
    }


        for(int i=0; i<m_nBin; i++) {
            double eTrue = hEdep->GetBinCenter(i+1);
            double tmp_entries = hEdep->GetBinContent(i+1);
            // nonlinearity
            double tmp_E = (electronQuench::ScintillatorPE(eTrue) + electronCerenkov::getCerPE(eTrue))/m_energyScale;
            // resolution
            //double tmp_sigmaE = electronResponse::fEvisSigma->Eval(tmp_E*m_energyScale) / m_energyScale;
            //double tmp_sigmaE = electronResponse::fEvisNew1->Eval(tmp_E*m_energyScale) / m_energyScale;
            double tmp_sigmaE = electronResponse::fEvisNew->Eval(tmp_E*m_energyScale) / m_energyScale;

            gaus->SetParameter(0, tmp_E);
            gaus->SetParameter(1, tmp_sigmaE);
            int minBin = int((tmp_E-5.0*tmp_sigmaE-m_EvisMin)/m_binCalcWidth);
            int maxBin = int((tmp_E+5.0*tmp_sigmaE-m_EvisMin)/m_binCalcWidth);

            for(int j=minBin; j<maxBin+1; j++) {
                if(j<0 or j>m_nBin) continue;
                double tmp_center = m_binCalcWidth /2 + m_binCalcWidth * j + m_EvisMin;
                double prob = gaus->Eval(tmp_center) * m_binCalcWidth;

                hTmpEvisCalc->AddBinContent(j+1, prob*tmp_entries);
            }
        }


        // Rebinning
        m_calcEntries = 0;
        for(int i=0; i<m_nBinData; i++) {
            double bin1    = hEvisCalc->GetBinCenter(i+1);
            double bin1low = bin1 - hEvisCalc->GetBinWidth(i+1)/2.;
            double bin1hig = bin1 + hEvisCalc->GetBinWidth(i+1)/2.;
            for(int j=0; j<m_nBin; j++) {
                double bin2    = hTmpEvisCalc->GetBinCenter(j+1);
                double entry   = hTmpEvisCalc->GetBinContent(j+1);
                if(bin2 > bin1hig)
                    break;
                if(bin2 > bin1low and bin2 <= bin1hig) {
                    hEvisCalc->AddBinContent(i+1, entry); 
                    m_calcEntries += entry;
                }
            }
        }


}

double michelElectron::GetChi2()
{
    CalcVisibleEnergy();

    double chi2 = 0;
    for(int i=0; i<m_nBinData; i++) {
        double E = hEvisData->GetBinError(i+1);
        if (E==0) continue;
        double P = hEvisCalc->GetBinContent(i+1) / m_calcEntries * m_dataEntries;
        double M = hEvisData->GetBinContent(i+1);
        chi2 += (P-M) * (P-M) / E / E;
    }

    return chi2;
}

void michelElectron::Plot()
{
    hEdep->SetLineColor(kBlue+1);
    hEdep->SetLineWidth(2);
    hEvisData->SetLineColor(kRed+1);
    hEvisData->SetLineWidth(2);
    hTmpEvisCalc->SetLineColor(kGreen+1);
    hTmpEvisCalc->SetLineWidth(2);
    hEvisCalc->SetLineColor(kOrange+1);
    hEvisCalc->SetLineWidth(2);

    hEvisCalc->Scale(m_dataEntries/m_calcEntries);

    TFile* out = new TFile("michel.root", "recreate");
    hEdep->Write();
    hEvisData->Write();
    hTmpEvisCalc->Write();
    hEvisCalc->Write();
    out->Close();
}
