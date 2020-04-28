#include "electronNLExperiment.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>

using namespace std;

bool electronNLExperiment::m_LoadData = false;
bool electronNLExperiment::m_CalcTheo = false;
double electronNLExperiment::m_energyScale = 1481.06;

TGraphErrors* electronNLExperiment::mTrueElectronNL;
TGraph*       electronNLExperiment::mFitElectronNL;

electronQuench* electronNLExperiment::mQuench = (electronQuench*)0;
electronCerenkov* electronNLExperiment::mCerenkov = (electronCerenkov*)0;

//electronNLExperiment::electronNLExperiment(electronQuench* aQuench, electronCerenkov* aCerenkov)
//{
//    mQuench = aQuench;
//    mCerenkov = aCerenkov;
//}

electronNLExperiment::electronNLExperiment()
{ /*cout << " >>> Initializing Naked Electron Data <<< " << endl;*/ }

electronNLExperiment::~electronNLExperiment()
{;}



void electronNLExperiment::LoadData ()  {
    cout << " >>> Loading Naked Electron Data <<< " << endl;
    mTrueElectronNL = new TGraphErrors();

    ifstream in;
    in.open(junoParameters::electronLSNL_File.c_str());
    if(!in){
        cout << " >>> Fail to Open Electron totPE File !! <<< " << endl;
    }
    string line;

    int count = 0;
    double tmp_E, tmp_totpe, tmp_sigma;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_E >> tmp_totpe >> tmp_sigma ;
        mTrueElectronNL->SetPoint(count, tmp_E, tmp_totpe/m_energyScale/tmp_E);
        mTrueElectronNL->SetPointError(count, 0, tmp_sigma*tmp_totpe/m_energyScale/tmp_E);
        count++;
    }
    in.close();

    m_LoadData = true;
}


void electronNLExperiment::UpdateDataElectronNL()
{
    if(!m_LoadData)  LoadData();
}


void electronNLExperiment::UpdateTheoElectronNL()
{
    if( !m_LoadData )  UpdateDataElectronNL();
    mFitElectronNL = new TGraph();

    const int data_num = mTrueElectronNL->GetN();
    double *Edep = mTrueElectronNL->GetX(); 
    double nl_pred[data_num];

    for(int i=0; i<data_num; i++){
        double fq = mQuench->ScintillatorNL(Edep[i]);
        double fc = mCerenkov->getCerenkovPE(Edep[i]);
        nl_pred[i] = fq + fc ;
        mFitElectronNL->SetPoint(i, Edep[i], nl_pred[i]);
        //cout << "TheoElecNL : " << Edep[i] << " " << fq << " " << fc << endl;
    }

    m_CalcTheo = true;
    
}

double electronNLExperiment::GetChi2 ( int nDoF )  {
    // every time should update theoNL
    m_CalcTheo = false;

    if (!m_LoadData) UpdateDataElectronNL();
    if (!m_CalcTheo) UpdateTheoElectronNL();
    double chi2 = 0;
    if( mTrueElectronNL->GetN() != mFitElectronNL->GetN() ){
        cout << " >>> Two Data have Different Lengths <<< " << endl; return -1.;
    } else {
        const int nPoints = mTrueElectronNL->GetN();
        double *dataNL = mTrueElectronNL->GetY();
        double *theoNL = mFitElectronNL->GetY();
        double *EvisError = mTrueElectronNL->GetEY();
        for(int iPoint=0; iPoint<nPoints; iPoint++){
            double error2 = EvisError[iPoint] * EvisError[iPoint];
            double delta  = dataNL[iPoint] - theoNL[iPoint];
            //cout << "dataNL: " << dataNL[iPoint] << "  theoNL: " << theoNL[iPoint] << 
            //"  error2:  " << error2 << "  delta_chi2 : " << TMath::Power(delta, 2)/error2 << endl;
            chi2 += TMath::Power(delta,2)/error2;
        }
        if( nDoF > 0 )  {
            chi2 /= double (nPoints - nDoF);
        }
        return chi2;
    }
}


void electronNLExperiment::Plot () {

    UpdateTheoElectronNL();

    cout << " ---> Plot Electron Fitting Results " << endl;
    mTrueElectronNL  ->SetLineColor(kBlue+1);
    mFitElectronNL   ->SetLineColor(kGreen+1);
    mTrueElectronNL  ->SetMarkerColor(kBlue+1);
    mFitElectronNL   ->SetMarkerColor(kGreen+1);
    mTrueElectronNL  ->SetMarkerStyle(20);
    mFitElectronNL   ->SetMarkerStyle(21);
    mTrueElectronNL  ->SetMarkerSize(1.1);
    mFitElectronNL   ->SetMarkerSize(1.1);

    TMultiGraph* fitGraph = new TMultiGraph();
    fitGraph->Add(mTrueElectronNL, "PZ");
    fitGraph->Add(mFitElectronNL,  "pZ");
    TCanvas* tmpC1 = new TCanvas();
    tmpC1->cd();  tmpC1->SetGrid();
    fitGraph->SetTitle("Electron NL Fitting; Etrue/MeV; Evis/Etrue");
    fitGraph->Draw("AP");
    //TCanvas* tmpC2 = new TCanvas();
    //tmpC2->cd();
    //mTrueElectronNL->Draw("APL");
    //mFitElectronNL ->Draw("PL SAME");
	TLegend* leg=new TLegend();
	leg->SetBorderSize(0);
	leg->SetFillColor(-1);
    leg->AddEntry(mTrueElectronNL, "Electron Data", "PE");
    leg->AddEntry(mFitElectronNL,  "Electron Fit" , "PE");
	leg->SetTextSize(19);
	leg->Draw("SAME");

    
    TFile* elecFile = new TFile (junoParameters::electronOut_File.c_str(), "recreate");
    tmpC1->Write("ElectronFit");
    elecFile->Close();
    delete elecFile;

}



