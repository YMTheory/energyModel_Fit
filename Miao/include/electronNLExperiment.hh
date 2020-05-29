#ifndef electronNLExperiment_h
#define electronNLExperiment_h

#include <vector>
#include <string>
#include "TGraph.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "junoExperiment.hh"
#include "junoParameters.hh"

using namespace std;

class electronQuench;
class electronCerenkov;

class electronNLExperiment //: public junoExperiment
{
    public:
        electronNLExperiment();
        //electronNLExperiment(electronQuench* aQuench, electronCerenkov* aCerenkov);
        ~electronNLExperiment();

        static void    setEnergyScale (double val) {m_energyScale = val;}
        static double  getEnergyScale ()           {return m_energyScale;}

        static void LoadData   ();
        static double GetChi2  (int nDoF = 0 );

        static void Plot       ();

    public:

        static void UpdateDataElectronNL();
        static void UpdateTheoElectronNL();
        //
        static TGraphErrors* GetTrueElectronNL() { return mTrueElectronNL; }
        static TGraph* GetFitElectronNL() { return mFitElectronNL; }
    
    private:
        static bool m_LoadData;
        static bool m_LoadSingleData;
        static bool m_LoadB12;
        static bool m_CalcTheo;
        static double m_energyScale;   // prior energy scale ...

        static TGraphErrors* mTrueElectronNL;
        static TGraph* mFitElectronNL;

        static electronQuench* mQuench;
        static electronCerenkov* mCerenkov;

        static vector<double> Etrue;

        static TH1D* mTrueB12Spec;
        static TH1D* mFitB12Spec;
        static TH1D* mTempB12Spec;
        static std::vector<double> predB12;
        static double B12_scale;
        
};

#endif
