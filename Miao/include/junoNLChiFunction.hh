#ifndef junoNLChiFunction_h
#define junoNLChiFunction_h

#include "junoB12Data.hh"
#include "junoC11Data.hh"
#include "junoC10Data.hh"
#include <TMinuit.h>

class junoNLChiFunction 
{
    public:
        junoNLChiFunction();
        ~junoNLChiFunction();

        void LoadData                 ();
        double GetChiSquare           ( double maxChi2 = 100000 );
        static void SetParameters     ( double *par );
        static double GetChi2         ( double maxChi2 = 100000 );  
        
        bool GridSearch               ();
        bool ScanContour              ();      
    
        static void Plot              ();

    private:
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

        TMinuit* junoNLMinuit;

        static double m_chi2;
        static double m_chi2Min;
        static double m_chi2Gam;
        static double m_chi2B12;
        static double m_chi2C11;
        static double m_chi2C10;
        static int m_nParameter;
        static double m_bestFit[20];
        static double m_bestFitError[20];
        static bool m_DoFit;
        static bool m_gridSearch;
        
        static junoB12Data* m_b12Data;
        static junoC11Data* m_c11Data;
        static junoC10Data* m_c10Data;

        static double final_kA;
        static double final_kB;
        static double final_kC;
        static double m_kALimit;
        static double m_kBLimit;
        static double m_kCLimit;

} ;


#endif
