#ifndef _GAMMANLEXPERIMENT_H
#define _GAMMANLEXPERIMENT_H

#include <vector>
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;

class electronQuench;
class electronCerenkov;
class primaryElectronDistribution;
class gammaNLExperiment
{
    public:
        gammaNLExperiment(electronQuench* aQuench, electronCerenkov* aCerenkov, primaryElectronDistribution* aPED);
        ~gammaNLExperiment();

    public:
        void Calculate(double *apar);

        void CalculateTrueGammaNL();
        void CalculateFitGammaNL(double *apar);

        double CalculateGammaNL(double *apar, string name);
        
        TGraphErrors* GetTrueGammaNL() { return mTrueGammaNL; }
        TGraph* GetFitGammaNL() { return mFitGammaNL; }

    private:
        TGraphErrors* mTrueGammaNL;
        TGraph* mFitGammaNL;

        static electronQuench* mQuench;
        static electronCerenkov* mCerenkov;
        static primaryElectronDistribution* mPED;

        vector<double> Etrue;
        vector<string> source_name;
};

#endif
