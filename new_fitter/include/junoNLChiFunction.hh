#ifndef junoNLChiFunction_h
#define junoNLChiFunction_h

#include "gammaData.hh"
#include "junoSpectrum.hh"

#include <TMinuit.h>

class junoNLChiFunction
{
    public:
        junoNLChiFunction();
        ~junoNLChiFunction();

        void LoadData();
        double GetChiSquare           ( double maxChi2 = 100000 );
        static void SetParameters     ( double *par );
        static double GetChi2         ( double maxChi2 = 100000 );  
    
        static void GammaPlot         ();
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
        

    private:
        static gammaData* Cs137Data;
        static gammaData* Mn54Data;
        static gammaData* nHData;
        static gammaData* K40Data;
        static gammaData* Co60Data;
        static gammaData* Tl208Data;
        static gammaData* nC12Data;
        static gammaData* O16Data;
        static gammaData* nFe56Data;

        static junoSpectrum* junoB12;

        static int m_nData;
        static std::string source_name[20];

        static gammaData* gammaData_array[20];

        static bool m_doGamFit;
        static bool m_doB12Fit;
};

#endif
