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
#include <TH1.h>

using namespace std;

bool electronNLExperiment::m_LoadData = false;
bool electronNLExperiment::m_LoadSingleData = false;
bool electronNLExperiment::m_LoadB12 = true;
bool electronNLExperiment::m_CalcTheo = false;
double electronNLExperiment::m_energyScale = 1481.06;

TGraphErrors* electronNLExperiment::mTrueElectronNL;
TGraph*       electronNLExperiment::mFitElectronNL;

TH1D* electronNLExperiment::mTrueB12Spec;
TH1D* electronNLExperiment::mFitB12Spec;
std::vector<double> electronNLExperiment::predB12;
double electronNLExperiment::B12_scale;

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

    if(m_LoadSingleData) {
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

    }


    if(m_LoadB12) {
        cout << " >>> Loading B12 Spectrum <<< " << endl;
        TFile* B12file = TFile::Open(junoParameters::B12DataFile.c_str());
        if(!B12file) cout << " >>> No B12 NL file <<< " << std::endl;
        mTrueB12Spec = (TH1D*)B12file->Get("Evis");
        mFitB12Spec  = (TH1D*)B12file->Get("Edep");
        
        ifstream in;
        in.open(junoParameters::B12PredFile.c_str());
        string line; double tmp_E;
        while(getline(in,line)){
            istringstream ss(line);
            ss >> tmp_E;
            predB12.push_back(tmp_E); 
        }

    }

    m_LoadData = true;
}


void electronNLExperiment::UpdateDataElectronNL()
{
    if(!m_LoadData)  LoadData();
}


void electronNLExperiment::UpdateTheoElectronNL()
{
    if( !m_LoadData )  UpdateDataElectronNL();

    if( m_LoadSingleData ) {  // single electron data (mock)
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
    }

    if( m_LoadB12 ) {  // B12 spectrum
        mFitB12Spec->Reset();
        for (int k=0; k<predB12.size();k++){
            double fq = mQuench->ScintillatorNL(predB12[k]);
            double fC = mQuench->ScintillatorNL(predB12[k]);
            double evis = fq*predB12[k];
            //double evis = (fq+fC)*predB12[k];
            mFitB12Spec->Fill(evis);
        }
        B12_scale = mTrueB12Spec->GetEntries() / mFitB12Spec->GetEntries();
        mFitB12Spec->Scale(B12_scale);
    }


    m_CalcTheo = true;
    
}

double electronNLExperiment::GetChi2 ( int nDoF )  {
    // every time should update theoNL
    m_CalcTheo = false;

    if (!m_LoadData) UpdateDataElectronNL();
    if (!m_CalcTheo) UpdateTheoElectronNL();
    double chi2 = 0;
    double nData = 0;
    
    if( m_LoadSingleData ) {
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
                nData++;
            }
        }
    }

    if( m_LoadB12 ) {
        const int nBins = mTrueB12Spec->GetNbinsX();
        for(int iBin=0; iBin<nBins; iBin++) {
            double delta = mTrueB12Spec->GetBinContent(iBin)-mFitB12Spec->GetBinContent(iBin);
            double error2 = mTrueB12Spec->GetBinCenter(iBin);
            if(error2!=0) {chi2 += TMath::Power(delta, 2)/error2; nData++;}
        }
    }
        if( nDoF > 0 )  {
            chi2 /= double (nData - nDoF);
        }
        return chi2;
}


void electronNLExperiment::Plot () {

    UpdateTheoElectronNL();

    TFile* elecFile = new TFile (junoParameters::electronOut_File.c_str(), "recreate");
    cout << " ---> Plot Electron Fitting Results " << endl;

    if( m_LoadSingleData ) {
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

        tmpC1->Write("ElectronFit");
    }

    if( m_LoadB12 ) {
        mTrueB12Spec  ->SetLineColor(kBlue+1);
        mFitB12Spec   ->SetLineColor(kGreen+1);
        mTrueB12Spec  ->SetLineWidth(2);
        mFitB12Spec   ->SetLineWidth(2);
        TCanvas* tmpC2 = new TCanvas();
        tmpC2->cd();  tmpC2->SetGrid();
        mTrueB12Spec->Draw();
        mFitB12Spec->Draw("SAMES");
        TLegend* leg=new TLegend();
        leg->SetBorderSize(0);
        leg->SetFillColor(-1);
        leg->AddEntry(mTrueB12Spec, "B12Spec Data", "PE");
        leg->AddEntry(mFitB12Spec,  "B12Spec Fit" , "PE");
        leg->SetTextSize(19);
        leg->Draw("SAME");
        
        tmpC2->Write("ElectronFit");
        
    }
    
    elecFile->Close();
    delete elecFile;

}



