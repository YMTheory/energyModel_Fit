#include "junoC10Data.hh"
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

unsigned int junoC10Data::s_nMaxBins     = 4000;
//unsigned int junoC10Data::s_nMaxBins     = 28800;
//unsigned int junoC10Data::s_nMaxBinsData = 80;
unsigned int junoC10Data::s_nMaxBinsData = 60;
unsigned int junoC10Data::s_nMaxBr       = 5;
unsigned int junoC10Data::s_nMaxGam      = 2;

void junoC10Data::SetParameters()
{
    m_eMinDraw        = 0.;
    m_eMin            = 0.;
    m_eMax            = 6.0;
    m_fitMin          = junoParameters::c10FitMinE;
    m_fitMax          = junoParameters::c10FitMaxE;
    m_name            = "C10";
    m_title           = "C10";
    is_positron       = true;
}

void junoC10Data::InitTheo  ()
{
    std::cout << " calculating theoretical C10 shape " << endl;
    TheoHistTree("data/electron/theo/C10_theo.root");
}

void junoC10Data::InitData (string fileName) 
{
	std::cout << " ----> Reading C10 data from " << fileName << std::endl;

	TFile* file = new TFile(fileName.c_str());
	//TH1F* sigH  = (TH1F*)file->Get("spec");
	TH1F* sigH  = (TH1F*)file->Get("c10");
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


void junoC10Data::Plot()
{
    cout << " >>> Plot C10 Fitting Results <<< " << endl;
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
    hData->GetXaxis()->SetRangeUser(m_fitMin, m_fitMax);
    hData->Draw("PEX0");
    hTheo->Draw("PEXO SAME");
    TLegend* ll = new TLegend();
    ll->AddEntry(hData, "B12 Data(Sim)", "l");
    ll->AddEntry(hTheo, "B12 Theo", "l");
    ll->Draw("SAME");
    tmpC->cd(2);
    gPad->SetTopMargin(small);
    gPad->SetTickx();
    hRela->GetXaxis()->SetRangeUser(m_fitMin, m_fitMax);
    hRela->Draw("PEX0");

    TFile* file = new TFile(junoParameters::C10Out_File.c_str(), "recreate");
    tmpC->Write();
    file->Close();

}


