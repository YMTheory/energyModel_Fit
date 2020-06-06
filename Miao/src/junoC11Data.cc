#include "junoC11Data.hh"
#include "junoParameters.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TPad.h>

using namespace std;

unsigned int junoC11Data::s_nMaxBins     = 4000;
//unsigned int junoC11Data::s_nMaxBins     = 28800;
//unsigned int junoC11Data::s_nMaxBinsData = 80;
unsigned int junoC11Data::s_nMaxBinsData = 40;
unsigned int junoC11Data::s_nMaxBr       = 5;
unsigned int junoC11Data::s_nMaxGam      = 2;

void junoC11Data::SetParameters()
{
	m_eMinDraw     = 0;
	m_eMin         = 0;
	m_eMax         = 3.0;
	m_fitMin       = junoParameters::c11FitMinE;
	m_fitMax       = junoParameters::c11FitMaxE;
	m_name         = "C11";
	m_title        = "C11";
    is_positron    = true;

}


void junoC11Data::InitTheo ()
{
	std::cout << " calculating theoretical C11 shape " << std::endl;
	TheoHistTree("data/electron/theo/C11_theo.root");
}


void junoC11Data::InitData(string fileName)
{
	std::cout << " ----> Reading C11 data from " << fileName << std::endl;

	TFile* file = new TFile(fileName.c_str());
	//TH1F* sigH  = (TH1F*)file->Get("spec");
	TH1F* sigH  = (TH1F*)file->Get("c11");
	for (int i=0;i!=m_nBinsData;i++)
	{
		double content = sigH->GetBinContent(i+1);
		double error   = sigH->GetBinError  (i+1);
		m_eData   [i] = content;
		m_eDataErr[i] = error;
	}
	delete sigH;
	file->Close();
	delete file;
    
}


void junoC11Data::Plot()
{
    cout << " >>> Plot C11 Fitting Results <<< " << endl;
    TH1D* hData = new TH1D("hData", "", m_nBinsData, m_eMin, m_eMax);
    TH1D* hTheo = new TH1D("hTheo", "", m_nBinsData, m_eMin, m_eMax);
    TH1D* hRela = new TH1D("hRela", "", m_nBinsData, m_eMin, m_eMax);

    hData->SetStats(0);
    hData->SetLineColor(kBlue+1);
    hData->SetLineWidth(2);
    hData->SetMarkerSize(0.8);
    hData->SetMarkerStyle(20);
    hData->SetMarkerColor(kBlue+1);
    hTheo->SetStats(0);
    hTheo->SetLineColor(kRed+1);
    hTheo->SetLineWidth(2);
    hTheo->SetMarkerSize(0.8);
    hTheo->SetMarkerStyle(20);
    hTheo->SetMarkerColor(kRed+1);
    hRela->SetStats(0);
    hRela->SetLineColor(kPink+2);
    hRela->SetLineWidth(2);
    hRela->SetMarkerSize(0.8);
    hRela->SetMarkerStyle(20);
    hRela->SetMarkerColor(kPink+2);


    for(int i=0; i<m_nBinsData; i++) {
        hData->SetBinContent(i+1, m_eData[i]);
        hTheo->SetBinContent(i+1, m_eTheo[i]);
        if(m_eTheo[i]!=0) {
            hRela->SetBinContent(i+1, m_eData[i]/m_eTheo[i]);
            hRela->SetBinError(i+1, m_eDataErr[i]/m_eTheo[i]);
        } else {
            hRela->SetBinContent(i+1, 0);
            hRela->SetBinError(i+1, 0);
        }
    }

    TCanvas * tmpC = new TCanvas();
    Float_t small = 1e-5;
    tmpC->Divide(1,2, small, small);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    tmpC->cd(1);
    gPad->SetBottomMargin(small);
    hData->GetYaxis()->SetRangeUser(0,12500);
    hData->Draw("PEX0");
    hTheo->Draw("PEX0 SAME");
    TLegend* ll = new TLegend();
    ll->AddEntry(hData, "C11 Data(Sim)", "l");
    ll->AddEntry(hTheo, "C11 Theo", "l");
    ll->Draw("SAME");
    tmpC->cd(2);
    gPad->SetTopMargin(small);
    gPad->SetTickx();
    hRela->Draw("PEX0");

    TFile* file = new TFile(junoParameters::C11Out_File.c_str(), "recreate");
    tmpC->Write();
    file->Close();

}

