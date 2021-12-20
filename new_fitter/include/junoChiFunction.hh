#ifndef junoChiFunction_h
#define junoChiFunction_h

#include "gammaResponse.hh"
#include "junoB12_simplified.hh"
#include "michelElectron.hh"

#include <TMinuit.h>

class junoChiFunction
{
    public:
        junoChiFunction();
        ~junoChiFunction();

        void LoadData();
        double GetChiSquare           ( double maxChi2 = 100000 );
        static void SetParameters     ( double *par );
        static double GetChi2         ( double maxChi2 = 100000 );  
        static void Plot              ();

        static double GetChi2Michel   ();

        static void Chi2Scanner       ();


    private:
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

        TMinuit* junoMinuit;

        static double m_chi2;
        static double m_chi2Min;
        static double m_chi2Gam;
        static double m_chi2B12;
        static int m_nParameter;
        static double m_bestFit[20];
        static double m_bestFitError[20];
        static bool m_DoFit;
        
        static junoB12_simplified* b12data;

        static michelElectron* michel;

        static gammaResponse* gamma_array[10];
        static string gamma_name[10];
        static double gamma_amp[10];
        static double init_amp[10];
    
        static int m_nData;
        static bool m_doGamFit;
        static bool m_doB12Fit;
        static bool m_doMichel;

        static double m_micNPE;
        static double m_micSPE;
        static double m_micNPEerr;
        static double m_micSPEerr;
        static double m_resp1;
        static double m_resp2;
};

#endif
