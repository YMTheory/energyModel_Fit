#ifndef _GAMMANLCHIFUNCTION_H
#define _GAMMANLCHIFUNCTION_H

#include <TMinuit.h>

class gammaNLExperiment;

class gammaNLChiFunction
{
    public:
        gammaNLChiFunction(gammaNLExperiment* aGammaNLExperiment);
        ~gammaNLChiFunction();

    public:
        double GetChiSquare();

    private:
        static gammaNLExperiment* mGammaNLExperiment;

        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

        TMinuit* gammaNLMinuit;

        static double errGamma;
};

#endif
