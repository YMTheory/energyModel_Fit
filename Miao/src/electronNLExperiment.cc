#include "electronNLExperiment.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>

using namespace std;

electronQuench* electronNLExperiment::mQuench = (electronQuench*)0;
electronCerenkov* electronNLExperiment::mCerenkov = (electronCerenkov*)0;

electronNLExperiment::electronNLExperiment(electronQuench* aQuench, electronCerenkov* aCerenkov)
{
    mQuench = aQuench;
    mCerenkov = aCerenkov;
}

electronNLExperiment::~electronNLExperiment()
{;}

void electronNLExperiment::Calculate(double *apar)
{
    CalculateTrueElectronNL();

    CalculateFitElectronNL(apar);
}

void electronNLExperiment::CalculateTrueElectronNL()
{
    mTrueElectronNL = new TGraphErrors();

    double A = 1481.06; // if we have a prior energy scale here ...

    ifstream in;
    in.open("/Users/yumiao/Documents/Works/github/energyModel_juno/nonlinearity/fit_param/data/electron_total.txt");
    if(!in){
        cout << " >>> Fail to Open Electron totPE File !! <<< " << endl;
    }
    string line;

    int count = 0;
    double tmp_E, tmp_totpe, tmp_sigma;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_E >> tmp_totpe >> tmp_sigma ;
        mTrueElectronNL->SetPoint(count, tmp_E, tmp_totpe/A/tmp_E);
        mTrueElectronNL->SetPointError(count, 0, tmp_sigma);
        count++;
    }
    in.close();
}

void electronNLExperiment::CalculateFitElectronNL(double *apar)
{
    mFitElectronNL = new TGraph();

    // model parameters
    double kA = apar[0];
    double kB = apar[1];
    double kC = apar[2];
    double A = 1481.06; // prior energy scale 

    const int data_num = mTrueElectronNL->GetN();
    double *Edep = mTrueElectronNL->GetX(); 
    double nl_pred[data_num];

    for(int i=0; i<data_num; i++){
        double fq = mQuench->Integral_BirkLaw(kB, Edep[i]);
        double fc = mCerenkov->getCerenkovPE(Edep[i]);
        nl_pred[i] = kA*fq + kC*fc/A/Edep[i];
        mFitElectronNL->SetPoint(i, Edep[i], nl_pred[i]);
    }
    
}

