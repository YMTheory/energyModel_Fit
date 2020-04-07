#ifndef junoGlobalFit_h
#define junoGlobalFit_h

#include "TMinuit.h"

#include "junoEnergyModel.hh"
#include "junoParameters.hh"
#include "junoGammaData.hh"

using namespace std;

class junoParameters;

class junoGlobalFit
{
    public:
        junoGlobalFit();
        ~junoGlobalFit();
        void          LoadData ();
        void          Fit      ();

        static void SetParameters  ();
        static double GetChi2      (double maxChi2=-1);
        
        static junoGammaData *m_gammaData;

    private:

        TMinuit*    m_minuit;
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);
    
        static double m_parameters[4];
        static int    m_nParameter;
        static int    m_nFreeParameter;

        static double m_chi2;
        static double m_chi2Gamma;

        static double m_chi2Min;

        static double m_bestFit     [20];
        static double m_bestFitError[20];
        static double m_covMatrix   [20][20];

};
#endif
