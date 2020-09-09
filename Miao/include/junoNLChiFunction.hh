#ifndef junoNLChiFunction_h
#define junoNLChiFunction_h

#include "gammaResol.hh"
#include "junoB12Data.hh"
#include "junoC11Data.hh"
#include "junoC10Data.hh"
#include <TMinuit.h>
#include <map>

class junoNLChiFunction 
{
    public:
        junoNLChiFunction();
        ~junoNLChiFunction();

        void LoadData                 ();
        double GetChiSquare           ( double maxChi2 = 100000 );
        static void SetParameters     ( double *par );
        static double GetChi2         ( double maxChi2 = 100000 );  
    
        static void Plot              ();
        static void GammaPlot         ();

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
        
        static gammaResol* Cs137Data;
        static gammaResol* Mn54Data;
        static gammaResol* nHData;
        static gammaResol* K40Data;
        static gammaResol* Co60Data;
        static gammaResol* Tl208Data;
        static gammaResol* nC12Data;
        static gammaResol* O16Data;
        static gammaResol* nFe56Data;

        static junoB12Data* m_b12Data;
        static junoC11Data* m_c11Data;
        static junoC10Data* m_c10Data;

        static double final_kA;
        static double final_kB;
        static double final_kC;
        static double m_kALimit;
        static double m_kBLimit;
        static double m_kCLimit;

        static int m_nData;
        static std::map<std::string, gammaResol*> mapGammaResol;
        static std::string source_name[20];

} ;


#endif
