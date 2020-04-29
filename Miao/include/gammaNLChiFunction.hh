#ifndef _GAMMANLCHIFUNCTION_H
#define _GAMMANLCHIFUNCTION_H

#include <TMinuit.h>

class gammaNLExperiment;

class gammaNLChiFunction
{
    public:
        gammaNLChiFunction();
        ~gammaNLChiFunction();

    public:
        double GetChiSquare         (double maxChi2 = 100000);
        static void SetParameters   (double *par);
        static double GetChi2       (double maxChi2 = 100000);

        static void   Plot           ();    
        static void DrawContour      (unsigned int N1, unsigned int N2);

    private:
        static gammaNLExperiment* mGammaNLExperiment;

        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

        TMinuit* gammaNLMinuit;

        static double errGamma;
        static double m_chi2;

        static int m_nParameter;
        static double m_bestFit[20];
        static double m_bestFitError[20];

        static bool m_DoFit;
};

#endif
