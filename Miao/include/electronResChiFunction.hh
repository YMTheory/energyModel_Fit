#ifndef _ELECTRONRESCHIFUNCTION_H
#define _ELECTRONRESCHIFUNCTION_H

#include "TMinuit.h"

class electronResChiFunction
{
    public:
        electronResChiFunction();
        ~electronResChiFunction();

    public:
        double GetChiSquare       (double maxChi2 = 1000000);
        static void SetParameters (double *par);
        static double GetChi2     (double maxChi2 = 1000000);

        static void Plot          ();

    private:
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);
        TMinuit* electronResMinuit;
        
        static double m_chi2;
        static double m_chi2Min;

        static int m_nParameter;
        static double m_bestFit[20];
        static double m_bestFitError[20];

        static bool m_DoFit;
};
#endif
