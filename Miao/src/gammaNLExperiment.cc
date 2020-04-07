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

using namespace std;

electronQuench* gammaNLExperiment::mQuench = (electronQuench*)0;
electronCerenkov* gammaNLExperiment::mCerenkov = (electronCerenkov*)0;
primaryElectronDistribution* gammaNLExperiment::mPED = (primaryElectronDistribution*)0;

gammaNLExperiment::gammaNLExperiment(electronQuench* aQuench, electronCerenkov* aCerenkov, primaryElectronDistribution* aPED)
{
    mQuench = aQuench;
    mCerenkov = aCerenkov;
    mPED = aPED;
}

gammaNLExperiment::~gammaNLExperiment()
{
    delete mQuench;
    delete mCerenkov;
    delete mPED;
}


void gammaNLExperiment::Calculate(double *apar)
{
    CalculateTrueGammaNL();

    CalculateFitGammaNL(apar);
}

void gammaNLExperiment::CalculateTrueGammaNL()
{
    mTrueGammaNL = new TGraphErrors();

    double A = 1481.06;  // prior electron energy scale ..

    ifstream in;
    in.open("./data/naked_gamma/singleGamma.txt");
    string line;

    int num = 0;
    string tmp_name; double tmp_E, tmp_totPE, tmp_totPE_err, tmp_sigma, tmp_sigma_err;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_name >> tmp_E >> tmp_totPE >> tmp_totPE_err >> tmp_sigma >> tmp_sigma_err;   ;
        source_name.push_back(tmp_name);
        mTrueGammaNL->SetPoint(num, tmp_E, tmp_totPE/A/tmp_E);
        mTrueGammaNL->SetPointError(num, 0, tmp_sigma/A/tmp_E);
        num++;
    }
    in.close();
}


void gammaNLExperiment::CalculateFitGammaNL(double *apar)
{
    mFitGammaNL = new TGraph();

    const int nPoint = mTrueGammaNL->GetN();
    double *Edep = mTrueGammaNL->GetX();
    double nl_pred[nPoint];

    for(int i=0; i<nPoint; i++){
        nl_pred[i] = CalculateGammaNL(apar, source_name[i]);
        mFitGammaNL->SetPoint(i, Edep[i], nl_pred[i]);
    }
}


double gammaNLExperiment::CalculateGammaNL(double *apar, string name)
{
    double A = 1481.06;

    // model parameters
    double kA = apar[0];
    double kB = apar[1];
    double kC = apar[2];

    
    mPED->read_distribution(name);
    const int nBins = mPED->getSize();
    if(nBins==0) return 0;
    double numerator = 0.; double denominator = 0.;
    for(int iBin=1; iBin<nBins-1; iBin++){

        double EE1 = mPED->getEtrue(iBin-1);
        double EE2 = mPED->getEtrue(iBin);

        double fq1 = mQuench->Integral_BirkLaw(kB, EE1) ;
        double fc1 = mCerenkov->getCerenkovPE(EE1);
        double electronNL1 = kA*fq1+kC*fc1/A/EE1;
        double fq2 = mQuench->Integral_BirkLaw(kB, EE2);
        double fc2 = mCerenkov->getCerenkovPE(EE2);
        double electronNL2 = kA*fq2+kC*fc2/A/EE2;

        numerator +=   (mPED->getCount(iBin-1)* EE1 *electronNL1 + mPED->getCount(iBin)*EE2*electronNL2 ) * (EE2-EE1) /2.;
        denominator += (mPED->getCount(iBin-1)*EE1 + mPED->getCount(iBin)*EE2)*(EE2-EE1)/2.;
    }

    if(denominator ==0) { cout << " >> Error Happens While CalculateGammaNL <<<" << endl; return 0;}
    return numerator/denominator;

}





