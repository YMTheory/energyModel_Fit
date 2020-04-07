#include <iostream>
#include <vector>

#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronNLExperiment.hh"
#include "electronNLChiFunction.hh"
#include "primaryElectronDistribution.hh"
#include "gammaNLExperiment.hh"
#include "gammaNLChiFunction.hh"

#include "TGraphErrors.h"
#include "TGraph.h"

using namespace std;

int main(){

    electronQuench* gelectronQuench = new electronQuench();
    //cout << "Check BirkLaw: " << gelectronQuench->Integral_BirkLaw(6.5e-3, 2) << endl;

    electronCerenkov* gelectronCerenkov = new electronCerenkov();

    electronNLExperiment* gelectronNLExperiment = new electronNLExperiment(gelectronQuench, gelectronCerenkov);

    electronNLChiFunction* chiFCN = new electronNLChiFunction(gelectronNLExperiment);
    double chi2 = chiFCN->GetChiSquare();
    cout << "chi2: " << chi2 << endl;

    //primaryElectronDistribution* gPED = new primaryElectronDistribution();

    //gammaNLExperiment* gGammaNLExperiment = new gammaNLExperiment(gelectronQuench, gelectronCerenkov, gPED);

    //gammaNLChiFunction* gGammeNLChiFCN = new gammaNLChiFunction(gGammaNLExperiment);
    //double chi2 = gGammeNLChiFCN->GetChiSquare();
    //cout << "chi2: " << chi2 << endl;
}


