#ifndef _ELECTRONNLCHIFUNCTION_H
#define _ELECTRONNLCHIFUNCTION_H

#include <TMinuit.h>

class electronNLExperiment;

class electronNLChiFunction
{

    public:
        electronNLChiFunction(electronNLExperiment* aelectronNLExperiment);
        ~electronNLChiFunction();

    public:
        double GetChiSquare();

    private:
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

        static electronNLExperiment* melectronNLExperiment;

        TMinuit* electronNLMinuit;

        static double errElectron;
};

#endif
