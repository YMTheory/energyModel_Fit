#ifndef gammaResChiFunction_h
#define gammaResChiFunction_h

#include "gammaResol.hh"
#include <map>
#include <TMinuit.h>

class gammaResChiFunction
{
    public:
        gammaResChiFunction();
        ~gammaResChiFunction();

        double LoadDate               ();
        double GetChiSquare           ( double maxChi2 = 100000 );
        static void SetParameters     ( double *par );
        static double GetChi2         ( double maxChi2 = 100000 );  

        static bool GridSearch        ();

        static void Plot              ();

        static void Contour           ();


    private:
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);
        TMinuit* gammaResMinuit;

        static double m_chi2;
        static int m_nParameter;
        static double m_bestFit[20];
        static double m_bestFitError[20];

        static gammaResol* Cs137Data;
        static gammaResol* Mn54Data;
        static gammaResol* nHData;
        static gammaResol* K40Data;
        static gammaResol* Co60Data;
        static gammaResol* Tl208Data;
        static gammaResol* nC12Data;
        static gammaResol* O16Data;
        static gammaResol* nFe56Data;

        static bool m_DoFit;
        static bool m_gridSearch;
        static int m_nData;
        static double final_pA;
        static double final_pB;
        static double final_pC;
        static double final_kA;
        static double final_kB;
        static double final_kC;

        static std::map<std::string, gammaResol*> mapGammaResol;
        static std::string source_name[20];

};

#endif
