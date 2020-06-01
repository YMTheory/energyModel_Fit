#include "junoB12Data.hh"
#include "junoParameters.hh"
#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;


unsigned int junoB12Data::s_nMaxBins     = 14400;
//unsigned int junoB12Data::s_nMaxBins     = 28800;
//unsigned int junoB12Data::s_nMaxBinsData = 80;
unsigned int junoB12Data::s_nMaxBinsData = 100;
unsigned int junoB12Data::s_nMaxBr       = 5;
unsigned int junoB12Data::s_nMaxGam      = 2;

void junoB12Data::SetParameters()
{
	m_eMinDraw     = 0;
	m_eMin         = 0;
	m_eMax         = 15;
	m_fitMin       = junoParameters::b12FitMinE;
	m_fitMax       = junoParameters::b12FitMaxE;
	m_name         = "B12";
	m_title        = "B12";
}

void junoB12Data::InitTheo ()
{
	std::cout << " calculating theoretical B12 shape " << std::endl;
	TheoHistTree("data/electron/theo/B12_hist.root");
}

void junoB12Data::InitData(string fileName)
{
	std::cout << " ----> Reading B12 data from " << fileName << std::endl;

	TFile* file = new TFile(fileName.c_str());
	//TH1F* sigH  = (TH1F*)file->Get("spec");
	TH1F* sigH  = (TH1F*)file->Get("b12");
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

void junoB12Data::Plot()
{
    cout << " >>> Plot B12 Fitting Results <<< " << endl;

    TH1D* hData = new TH1D("hData", "", m_nBinsData, m_eMin, m_eMax);
    TH1D* hTheo = new TH1D("hTheo", "", m_nBinsData, m_eMin, m_eMax);

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


    for(int i=0; i<m_nBinsData; i++) {
        hData->SetBinContent(i+1, m_eData[i]);
        hTheo->SetBinContent(i+1, m_eTheo[i]);
    }

    TCanvas * tmpC = new TCanvas();
    tmpC->cd();
    hData->Draw("PEX0");
    hTheo->Draw("PEXO SAME");
    TLegend* ll = new TLegend();
    ll->AddEntry(hData, "B12 Data(Sim)", "l");
    ll->AddEntry(hTheo, "B12 Theo", "l");
    ll->Draw("SAME");
    TFile* file = new TFile(junoParameters::B12Out_File.c_str(), "recreate");
    tmpC->Write();
    file->Close();
}
