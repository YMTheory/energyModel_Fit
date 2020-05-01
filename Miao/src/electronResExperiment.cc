#include "electronResExperiment.hh"
#include "junoParameters.hh"
#include "electronResol.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

electronResExperiment::electronResExperiment()  {;}
electronResExperiment::~electronResExperiment() {B12file->Close(); C11file->Close();}

bool electronResExperiment::m_LoadData = false;
bool electronResExperiment::m_LoadB12Data = true;
bool electronResExperiment::m_LoadC11Data = true;
bool electronResExperiment::m_CalcTheo = false;

double electronResExperiment::m_B12Energy[500] = {0.};
double electronResExperiment::m_B12Origin[500] = {0.};
double electronResExperiment::m_B12Nonl[500]   = {0.};
double electronResExperiment::m_B12Smear[500]  = {0.};
std::vector<double> electronResExperiment::B12energySpec; 
std::vector<double> electronResExperiment::C11energySpec; 
TH1D* electronResExperiment::h_B12Origin;
TH1D* electronResExperiment::h_B12Nonl;
TH1D* electronResExperiment::h_B12Smear;
TH1D* electronResExperiment::h_B12Fit;

TFile* electronResExperiment::B12file;
TFile* electronResExperiment::C11file;

TH1D* electronResExperiment::h_C11Origin;
TH1D* electronResExperiment::h_C11Nonl;
TH1D* electronResExperiment::h_C11Smear;
TH1D* electronResExperiment::h_C11Fit;

void electronResExperiment::LoadData() {
    if(m_LoadB12Data) LoadB12Data();
    if(m_LoadC11Data) LoadC11Data();
    m_LoadData  = true; 
}

void electronResExperiment::LoadB12Data() {
    cout << " >>> Load B12 Data <<< " << endl;
    //TFile* file = TFile::Open(junoParameters::B12_File.c_str());
    B12file = TFile::Open(junoParameters::B12_File.c_str());
    h_B12Origin = (TH1D*)B12file->Get("B12Origin");
    h_B12Nonl   = (TH1D*)B12file->Get("B12Nonl");
    h_B12Smear  = (TH1D*)B12file->Get("B12Smear");
    h_B12Fit = (TH1D*)h_B12Nonl->Clone("B12Fit");
    
    //const int nBins = h_B12Nonl->GetNbinsX(); 
    //for(int iBin=0; iBin<nBins; iBin++) {
    //    m_B12Energy[iBin] = h_B12Origin->GetBinCenter(iBin);
    //    m_B12Origin[iBin] = h_B12Origin->GetBinContent(iBin);
    //    m_B12Nonl[iBin]   = h_B12Nonl->GetBinContent(iBin); 
    //    m_B12Smear[iBin]  = h_B12Smear->GetBinContent(iBin);
    //}

    // read energy spectrum:
    ifstream in;
    in.open(junoParameters::B12Spec_File.c_str());
    string line; double tmp; int num = 0;
    B12energySpec.clear();
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> tmp;  B12energySpec.push_back(tmp);
    }
    in.close();

}

void electronResExperiment::LoadC11Data() {
    cout << " >>> Load C11 Data <<< " << endl;

    C11file = TFile::Open(junoParameters::C11_File.c_str());
    h_C11Origin = (TH1D*)C11file->Get("C11Origin");
    h_C11Nonl   = (TH1D*)C11file->Get("C11Nonl");
    h_C11Smear  = (TH1D*)C11file->Get("C11Smear");
    h_C11Fit = (TH1D*)h_C11Nonl->Clone("C11Fit");
    

    // read energy spectrum:
    ifstream in;
    in.open(junoParameters::C11Spec_File.c_str());
    string line; double tmp; int num = 0;
    C11energySpec.clear();
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> tmp;  C11energySpec.push_back(tmp);
    }
    in.close();

}

void electronResExperiment::UpdateTrueSpectrum() {
    if(!m_LoadData) LoadData();
}

void electronResExperiment::UpdateTheoSpectrum() {
    UpdateTrueSpectrum();
    
    if(m_LoadB12Data) {
        h_B12Fit->Reset();
        const int vecSizeB12 = B12energySpec.size();
        for(int i=0; i<vecSizeB12; i++) {
            double Esmear = electronResol::energySmear(B12energySpec[i]);
            h_B12Fit->Fill(Esmear);
        }
    }

    if(m_LoadC11Data) {
        h_C11Fit->Reset();

        const int vecSizeC11 = C11energySpec.size();
        for(int i=0; i<vecSizeC11; i++) {
            double Esmear = electronResol::energySmear(C11energySpec[i]);
            h_C11Fit->Fill(Esmear);
        }

    }

    m_CalcTheo = true;
}

double electronResExperiment::GetChi2() {
    m_CalcTheo = false;
    if(!m_CalcTheo) UpdateTheoSpectrum();
    double chi2 = 0;
    if(m_LoadB12Data) { 
        if(h_B12Smear->GetNbinsX() != h_B12Fit->GetNbinsX() ) {
            cout << " >>> B12 Fitting and Data have different binning !!! <<< " << endl; return -1;
        }
        const int nBins = h_B12Smear->GetNbinsX();
        for(int iBin =0; iBin<nBins; iBin++) {
            double delta = h_B12Smear->GetBinContent(iBin) - h_B12Fit->GetBinContent(iBin);
            double sigma = h_B12Smear->GetBinError(iBin);
            if(sigma!=0) chi2+=delta*delta/sigma/sigma; ;
        } 
    }

    if(m_LoadC11Data) { 
        if(h_C11Smear->GetNbinsX() != h_C11Fit->GetNbinsX() ) {
            cout << " >>> c11 Fitting and Data have different binning !!! <<< " << endl; return -1;
        }
        const int nBins = h_C11Smear->GetNbinsX();
        for(int iBin =0; iBin<nBins; iBin++) {
            double delta = h_C11Smear->GetBinContent(iBin) - h_C11Fit->GetBinContent(iBin);
            double sigma = h_C11Smear->GetBinError(iBin);
            if(sigma!=0) chi2+=delta*delta/sigma/sigma;
        }
    }

        return chi2;
}


void electronResExperiment::PlotB12() {
    if(!m_LoadB12Data) return;
    cout << " >>> Plot B12 Spectrum <<< " << endl;
    if(!m_CalcTheo) UpdateTheoSpectrum();

    TCanvas* cc = new TCanvas();
    h_B12Origin->SetLineColor(kRed+1);
    h_B12Nonl->SetLineColor(kBlue+1);
    h_B12Smear->SetLineColor(kGreen+1);
    h_B12Fit->SetLineColor(kViolet);

    h_B12Origin->SetTitle("B12 Beta Spectrum; energy/MeV; ");
    h_B12Origin->Draw();
    h_B12Nonl->Draw("SAME");
    h_B12Smear->Draw("SAME");
    h_B12Fit->Draw("SAME");
    h_B12Fit->Draw("SAME");

    TLegend* led = new TLegend();
    led->AddEntry(h_B12Origin, "origin", "l");
    led->AddEntry(h_B12Nonl, "nonlinearity", "l");
    led->AddEntry(h_B12Smear, "smear", "l");
    led->AddEntry(h_B12Fit, "fitting", "l");
    led->SetFillColor(kWhite);
    led->Draw("SAME");
    
    TFile* file = new TFile(junoParameters::B12Out_File.c_str(), "RECREATE");
    cc->Write("B12Fitting");
    file->Close();
}


void electronResExperiment::PlotC11() {
    if(!m_LoadC11Data) return;
    cout << " >>> Plot C11 Spectrum <<< " << endl;
    if(!m_CalcTheo) UpdateTheoSpectrum();

    TCanvas* cc = new TCanvas();
    h_C11Origin->SetLineColor(kRed+1);
    h_C11Nonl->SetLineColor(kBlue+1);
    h_C11Smear->SetLineColor(kGreen+1);
    h_C11Fit->SetLineColor(kViolet);

    h_C11Origin->SetTitle("C11 Beta+ Spectrum; energy/MeV; ");
    h_C11Origin->Draw();
    h_C11Nonl->Draw("SAME");
    h_C11Smear->Draw("SAME");
    h_C11Fit->Draw("SAME");

    TLegend* led = new TLegend();
    led->AddEntry(h_C11Origin, "origin", "l");
    led->AddEntry(h_C11Nonl, "nonlinearity", "l");
    led->AddEntry(h_C11Smear, "smear", "l");
    led->AddEntry(h_C11Fit, "fitting", "l");
    led->SetFillColor(kWhite);
    led->Draw("SAME");
    
    TFile* file = new TFile(junoParameters::C11Out_File.c_str(), "RECREATE");
    cc->Write("C11Fitting");
    file->Close();
}

