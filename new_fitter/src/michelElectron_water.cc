#include "michelElectron_water.hh"

#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TH1.h>
#include <string>

michelElectron_water::michelElectron_water()
{
    m_nBin = 7000;
    m_nBinData = 100;

    m_eMin = 0;
    m_eMax = 70;
    m_peMin = 0;
    m_peMax = 1540;
    m_peMinFit = 700;
    m_peMaxFit = 900;

    m_c0 = 15.6;
    m_p1 = 2.0;
    m_p2 = 0.001546;

    fCerNPE_water = new TF1("fCerNPE_water", "[0]*x", m_eMin, m_eMax);
    fCerPESigma_water = new TF1("fCerPESigma_water", "TMath::Sqrt([0]*x + [1]*x*x)", m_peMin, m_peMax);

    gCerNPE_water = new TGraph();
}

michelElectron_water::~michelElectron_water()
{
    delete fCerNPE_water;
    delete fCerPESigma_water;
    delete gCerNPE_water;
}

void michelElectron_water::Initialize()
{
    cout << " >>>>> Initializing water phase michel electrons spectrum <<<<< " << endl; 

    m_eBinWidth = (m_eMax - m_eMin) / m_nBin;
    m_peBinWidth = (m_peMax - m_peMin) / m_nBin;
    for(int i=0; i<m_nBin; i++) {
        m_binCenter[i] = m_eBinWidth/2. + m_eBinWidth * i;
    }
    cout << "Energy range binWidth = " << m_eBinWidth << ", NPE range binWidth = "<< m_peBinWidth << endl;

    LoadDataSpec();
    LoadTheoSpec();
    LoadSimNPE();

    gaus = new TF1("gaus", "1/(TMath::Sqrt(2*TMath::Pi())*[1]) * TMath::Exp(-(x-[0])*(x-[0])/2/[1]/[1])", -100, 30000);

}


void michelElectron_water::LoadDataSpec()
{
    TH1D* sigH = new TH1D("michel_data", "", m_nBinData, m_peMin, m_peMax);

    TFile* ff = new TFile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/water/michel/michel_totpe_water1.root", "read");
    TTree* tt = (TTree*)ff->Get("michel");
    int m_totpe;
    tt->SetBranchAddress("totpe", &m_totpe);

    for(int i=0; i<tt->GetEntries(); i++) {
        tt->GetEntry(i);
        sigH->Fill(m_totpe);
    }

    for (int i=0; i<m_nBinData; i++) {
		double content = sigH->GetBinContent(i+1);
		double error   = sigH->GetBinError  (i+1);
		m_eData   [i] = content;
		m_eDataErr[i] = error;
    }
    
	delete sigH;
    
    delete tt;
    delete ff;
    m_loadData = true;
    cout << " >>> Load Data Spectrum for water phase michel electrons <<< " << endl;

}

void michelElectron_water::LoadTheoSpec()
{

    TH1F* simH;

    TFile* ff = new TFile("./data/electron/michel/michel_theo.root", "read");
    simH = (TH1F*)ff->Get("hh0");

    for(int i=0; i<m_nBin; i++) {
        m_eTru[i] = simH->GetBinContent(i+1);
    }

    delete simH;
    delete ff;
    m_loadTheo = true;
    cout << " >>> Load Theo Spectrum for water phase michel electrons <<< " << endl;
}


void michelElectron_water::LoadSimNPE()
{
    ifstream in; in.open("./data/electron/cerPESigma_water.txt");
    string line;
    double E, npe, npe_sigma;
    int cnt = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> E >> npe >> npe_sigma;
        gCerNPE_water->SetPoint(cnt, E, npe);
        cnt++;
    }

    in.close();

}



void michelElectron_water::ApplyLSResponse()
{
    if (not m_loadData) LoadDataSpec();
    if (not m_loadTheo) LoadTheoSpec();

    SetParameters();

    int newBin, newBinLow, newBinHig;
    double bias;

    for (int i=0; i<m_nBin; i++) {
        m_eVis[i] = 0;
    }

    for (int i=0; i<m_nBin; i++) {
        double eTru = m_binCenter[i];
        //double tmp_pe = fCerNPE_water->Eval(eTru); 
        double tmp_pe = gCerNPE_water->Eval(eTru);
        double tmp_sigma = fCerPESigma_water->Eval(tmp_pe);

        gaus->SetParameter(0, tmp_pe);
        gaus->SetParameter(1, tmp_sigma);

        int minBin = int((tmp_pe-3*tmp_sigma)/m_peBinWidth);
        int maxBin = int((tmp_pe+3*tmp_sigma)/m_peBinWidth);
        //if (minBin < 0) minBin = 0;
        //if(maxBin >= m_nBin) maxBin = m_nBin-1;

        for (int j=minBin; j<maxBin; j++) {
            if (j<0 or j>=m_nBin) continue;
            double tmp_center = m_peBinWidth /2 + m_peBinWidth * j;
            double prob = gaus->Eval(tmp_center);

            m_eVis[j] += prob * m_eTru[i];
        }

    }

}

void michelElectron_water::Normalize()
{
    int rebin = m_nBin / m_nBinData;
    double binWidthData = m_peBinWidth * rebin;
    double nTheo = 0;
    double nData = 0;
    for (int i=0; i<m_nBinData; i++) {
        m_eTheo[i] = 0;
        for (int j=0; j<rebin; j++) {
            m_eTheo[i] += m_eVis[i*rebin+j];
        }
        if (i*binWidthData>m_peMinFit and i*binWidthData<m_peMaxFit) {
            nTheo += m_eTheo[i];
            nData += m_eData[i];
        }
    }

    double scale = 1;
    if( nTheo!=0 ) { scale = nData/nTheo; }
	for (int i = 0; i < m_nBinData; i++)
    {
		m_eTheo[i] *= scale;
        
    }
	for (int i = 0; i < m_nBinData; i++)
	{
		m_eVis   [i] *= scale;
	}

}


double michelElectron_water::GetChi2()
{
    ApplyLSResponse();
    Normalize();

    double chi2 = 0;
    int rebin = m_nBin / m_nBinData;
	double binWidthData = m_peBinWidth * rebin;
    int m_nData = 0;
    for(int i=0; i < m_nBinData; i++) {
        if(i*binWidthData<m_peMinFit or binWidthData*i>m_peMaxFit) continue;
        if( m_eDataErr[i]!=0 ) {
            //cout << m_eData[i] << " " << m_eTheo[i] << " " << m_eDataErr[i] << endl;
            chi2 += pow( (m_eData[i] - m_eTheo[i])/m_eDataErr[i], 2); 
            m_nData++;
        }
    }
    //cout << "simplified B12 chi2: " << chi2 << " with nData : " << m_nData << endl;
	//if(nDoF>0) chi2 /= double(m_nData - nDoF);
	return chi2;

}


void michelElectron_water::Plot()
{
        
    TH1D* hData = new TH1D("hData", "", m_nBinData, m_peMin, m_peMax);
    TH1D* hTheo = new TH1D("hTheo", "", m_nBinData, m_peMin, m_peMax);
    TH1D* hRela = new TH1D("hRela", "", m_nBinData, m_peMin, m_peMax);

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

    for(int i=0; i<100; i++) {
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


    TFile* out = new TFile("michel_spectrum.root", "recreate");
    hData->Write();
    hTheo->Write();
    hRela->Write();
    out->Close();

}


void michelElectron_water::SetParameters()
{
    fCerNPE_water->SetParameter(0, m_c0);

    fCerPESigma_water->SetParameter(0, m_p1);
    fCerPESigma_water->SetParameter(1, m_p2);
}



void michelElectron_water::Compare()
{
    // set a series of parameters input :
    //
    TH1D* hData = new TH1D("hData", "", m_nBinData, m_peMin, m_peMax);
    for(int i=0; i<100; i++) {
        hData->SetBinContent(i+1, m_eData[i]);
    }
    hData->SetStats(0);
    hData->SetLineColor(kBlue+1);
    hData->SetLineWidth(2);
    hData->SetMarkerSize(0.8);
    hData->SetMarkerStyle(20);
    hData->SetMarkerColor(kBlue+1);

    TH1D* hComp[9];
    string name[9] = {"comp1", "comp2", "comp3", "comp4", "comp5", "comp6", "comp7", "comp8", "comp9"};


    double tmp_p1 = 1.76105;
    double tmp_p2 = 0.00197229;
    int cnt = 0;
    for(int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            m_p1 = tmp_p1 * (0.9 + 0.1 * i); 
            m_p2 = tmp_p2 * (0.9 + 0.1 * j); 
            SetParameters();
            cout << "current chi2 = " << GetChi2() << " , for p1 = " << m_p1 << " p2 = " << m_p2 << endl;

            string histname = name[cnt];
            hComp[cnt] = new TH1D(histname.c_str(), "", m_nBinData, m_peMin, m_peMax);
            for(int k=0; k<100; k++) {
                hComp[cnt]->SetBinContent(k+1, m_eTheo[k]);
            }
            hComp[cnt]->SetLineColor(20+cnt*2);
            cnt++;
        }
    }

    
    TFile* out = new TFile("michel_resol_compare.root", "recreate");
    hData->Write();
    for(int i=0; i<9; i++) {
        hComp[i]->Write();
    }
    out->Close();



}









