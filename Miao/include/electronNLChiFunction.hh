#ifndef _ELECTRONNLCHIFUNCTION_H
#define _ELECTRONNLCHIFUNCTION_H

#include <TMinuit.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronNLExperiment.hh"

class electronNLExperiment;

class electronNLChiFunction
{

    public:
        electronNLChiFunction();
        ~electronNLChiFunction();

    public:
        double        GetChiSquare   (double maxChi2 = 10000);
        static void   SetParameters  (double *par);
        static double GetChi2        (double maxChi2 = 10000);

        static void   Plot           ();    

        static void DrawContour      (unsigned int N1, unsigned int N2);

    private:

        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

        TMinuit* electronNLMinuit;

        static double errElectron;
        static double m_chi2;
        static double m_chi2Min;

        static int m_nParameter;
        static double m_bestFit[20];
        static double m_bestFitError[20];

        static bool m_DoFit;
};

#endif
