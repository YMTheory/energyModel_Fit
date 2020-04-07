#ifndef electronNLExperiment_h
#define electronNLExperiment_h

#include <vector>
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;

class electronQuench;
class electronCerenkov;

class electronNLExperiment
{
    public:
        electronNLExperiment(electronQuench* aQuench, electronCerenkov* aCerenkov);
        ~electronNLExperiment();

    public:
        void Calculate(double *apar);

        void CalculateTrueElectronNL();
        void CalculateFitElectronNL(double *apar);
        
        TGraphErrors* GetTrueElectronNL() { return mTrueElectronNL; }
        TGraph* GetFitElectronNL() { return mFitElectronNL; }
    
    private:
        TGraphErrors* mTrueElectronNL;
        TGraph* mFitElectronNL;

        static electronQuench* mQuench;
        static electronCerenkov* mCerenkov;

        vector<double> Etrue;
        
};

#endif
