#include "gammaNLExperiment.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "primaryElectronDistribution.hh"

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

bool gammaNLExperiment::m_LoadGammaData = false;
bool gammaNLExperiment::m_LoadPrmElecDist = false;
bool gammaNLExperiment::m_CalcTheo = false;

double gammaNLExperiment::m_gammaScale = 0.0;

double gammaNLExperiment::m_pdf_eTru[m_nMaxSources][m_nMaxPdf] = {0.};
double gammaNLExperiment::m_pdf_prob[m_nMaxSources][m_nMaxPdf] = {0.};
double gammaNLExperiment::m_max_eTru[m_nMaxSources] = {0.};

vector<double> gammaNLExperiment::Etrue;
vector<string> gammaNLExperiment::source_name;

TGraphErrors* gammaNLExperiment::mTrueGammaNL;
TGraph*       gammaNLExperiment::mFitGammaNL;

//gammaNLExperiment::gammaNLExperiment(electronQuench* aQuench, electronCerenkov* aCerenkov, primaryElectronDistribution* aPED)
//{
//    mQuench = aQuench;
//    mCerenkov = aCerenkov;
//    mPED = aPED;
//}

gammaNLExperiment::gammaNLExperiment()
{;}

gammaNLExperiment::~gammaNLExperiment()
{;}


void gammaNLExperiment::Calculate(double *apar)
{
    //CalculateTrueGammaNL();

    //CalculateFitGammaNL(apar);
}

void gammaNLExperiment::LoadData () {
    cout << " >>> Loading Naked Gamma Data <<< " << endl;
    mTrueGammaNL = new TGraphErrors();

    ifstream in;
    in.open(junoParameters::gammaLSNL_File);
    string line;

    int num = 0;
    string tmp_name; double tmp_E, tmp_Evis, tmp_EvisError;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_name >> tmp_E >> tmp_Evis >> tmp_EvisError;
        source_name.push_back(tmp_name);
        //mTrueGammaNL->SetPoint(num, tmp_E, (1+m_gammaScale)*tmp_Evis/tmp_E);
        //mTrueGammaNL->SetPointError(num, 0, (1+m_gammaScale)*tmp_EvisError*tmp_Evis/tmp_E);
        mTrueGammaNL->SetPoint(num, tmp_E, tmp_Evis/tmp_E);
        mTrueGammaNL->SetPointError(num, 0, tmp_EvisError*tmp_Evis/tmp_E);
        num++;
    }
    in.close();

    m_LoadGammaData = true;
}


void gammaNLExperiment::UpdateDataGammaNL()
{
    if(!m_LoadGammaData) LoadData();
}

void gammaNLExperiment::LoadPrimaryElecDist ()  {

    if( !m_LoadGammaData )  UpdateDataGammaNL();  // load gamma data firstly ...

    TFile file(junoParameters::gammaPdf_File.c_str(), "read");
    const int nSources = source_name.size();
    for (int iSource = 0; iSource<nSources; iSource++) {
        string pdfName = "gamma"+source_name[iSource];
        TGraph* gGammaPdf = (TGraph*)file.Get(pdfName.c_str());
        if(!gGammaPdf) cout << "No Such Pdf : " << pdfName << endl;
        
        // make sure all sources with data have corresponding pdf ... 
        double *tmp_eTru = gGammaPdf->GetX();
        double *tmp_prob = gGammaPdf->GetY();
        for(int i=0; i<gGammaPdf->GetN(); i++)  {
            m_pdf_eTru[iSource][i] = tmp_eTru[i];
            m_pdf_prob[iSource][i] = tmp_prob[i];
            if (m_pdf_eTru[iSource][i] > 0) m_max_eTru[iSource] = i;    // max energy cut ... 
        }
    }
    file.Close();
    m_LoadPrmElecDist = true;  
}



void gammaNLExperiment::UpdateTheoGammaNL()
{
    if( !m_LoadPrmElecDist ) LoadPrimaryElecDist();
    mFitGammaNL = new TGraph();

    const int nSources = mTrueGammaNL->GetN();
    double *Edep = mTrueGammaNL->GetX();
    double nl_pred[nSources];

    for(int iSource=0; iSource<nSources; iSource++){
        nl_pred[iSource] = CalculateGammaNL(iSource);
        mFitGammaNL->SetPoint(iSource, Edep[iSource], nl_pred[iSource]);
    }

    m_CalcTheo = true;
}


double gammaNLExperiment::CalculateGammaNL( int iSource )
{
    if( !m_LoadPrmElecDist ) LoadPrimaryElecDist();
    
    double numerator = 0.; double denominator = 0.;
    for(int iBin=1;  iBin<m_max_eTru[iSource]; iBin++) {
        double E1 = m_pdf_eTru[iSource][iBin-1];
        double E2 = m_pdf_eTru[iSource][iBin];

        double prob1 = m_pdf_prob[iSource][iBin-1];
        double prob2 = m_pdf_prob[iSource][iBin];

        double fNL1 = electronQuench::ScintillatorNL(E1) + electronCerenkov::getCerenkovPE(E1);
        double fNL2 = electronQuench::ScintillatorNL(E2) + electronCerenkov::getCerenkovPE(E2);

        numerator   += ( prob1*E1*fNL1 + prob2*E2*fNL2 ) * (E2-E1) /2.;
        denominator += (prob1*E1 + prob2*E2) * (E2-E1)/ 2.;
    } 

    if(denominator ==0) { cout << " >> Error Happens While CalculateGammaNL <<<" << endl; return 0;}
    return numerator/denominator;

}


double gammaNLExperiment::GetChi2 ( int nDoF ) {
    m_CalcTheo = false;


    if (!m_LoadGammaData)      UpdateDataGammaNL();
    if (!m_LoadPrmElecDist)    LoadPrimaryElecDist();
    if (!m_CalcTheo)           UpdateTheoGammaNL();

    double chi2 = 0;

    if( mTrueGammaNL->GetN() != mFitGammaNL->GetN() ){
        cout << " >>> Two Data have Different Lengths <<< " << endl; return -1.;
    } else {
        const int nPoints = mTrueGammaNL->GetN();
        double *dataNL = mTrueGammaNL->GetY();
        double *theoNL = mFitGammaNL->GetY();
        double *EvisError = mTrueGammaNL->GetEY();
        for(int iPoint=0; iPoint<nPoints; iPoint++){
            double error2 = EvisError[iPoint] * EvisError[iPoint];
            double delta  = dataNL[iPoint] - theoNL[iPoint];
            //cout << "dataNL: " << dataNL[iPoint] << "  theoNL: " << theoNL[iPoint] << endl;
            //"  error2:  " << error2 << "  delta_chi2 : " << TMath::Power(delta, 2)/error2 << endl;
            chi2 += TMath::Power(delta,2)/error2;
        }
        if( nDoF > 0 )  {
            chi2 /= double (nPoints - nDoF);
        }
        //cout << "GetGammaChi2 : " << chi2 << endl;
        return chi2;
    }
}

void gammaNLExperiment::Plot ()  {
    UpdateTheoGammaNL ();

    cout << " ---> Plot Gamma Fitting Results " << endl;
    mTrueGammaNL  ->SetLineColor(kBlue+1);
    mFitGammaNL   ->SetLineColor(kGreen+1);
    mTrueGammaNL  ->SetMarkerColor(kBlue+1);
    mFitGammaNL   ->SetMarkerColor(kGreen+1);
    mTrueGammaNL  ->SetMarkerStyle(20);
    mFitGammaNL   ->SetMarkerStyle(21);
    mTrueGammaNL  ->SetMarkerSize(1.1);
    mFitGammaNL   ->SetMarkerSize(1.1);
    mTrueGammaNL  ->SetLineWidth(2);
    mFitGammaNL   ->SetLineWidth(2);

    TMultiGraph* fitGraph = new TMultiGraph();
    fitGraph->Add(mTrueGammaNL, "PZ");
    fitGraph->Add(mFitGammaNL,  "pZ");
    TCanvas* tmpC1 = new TCanvas();
    tmpC1->cd();  tmpC1->SetGrid();
    fitGraph->SetTitle("Gamma NL Fitting; Etrue/MeV; Evis/Etrue");
    fitGraph->Draw("APL");
    //TCanvas* tmpC2 = new TCanvas();
    //tmpC2->cd();
    //mTrueGammaNL->Draw("APL");
    //mFitGammaNL ->Draw("PL SAME");
	TLegend* leg=new TLegend();
	leg->SetBorderSize(0);
	leg->SetFillColor(-1);
    leg->AddEntry(mTrueGammaNL, "Gamma Data", "PE");
    leg->AddEntry(mFitGammaNL,  "Gamma Fit" , "PE");
	leg->SetTextSize(19);
	leg->Draw("SAME");

    
    TFile* gammaFile = new TFile (junoParameters::gammaOut_File.c_str(), "recreate");
    tmpC1->Write("GammaFit");
    gammaFile->Close();
    delete gammaFile;
}



