#include "global.hh"

#include <TF1.h>

#ifndef radiative_h
#define radiative_h

class radiative()
{
    public:
        radiative();
        double radiaCorr(double betaE);
        void Plot();
    private:
        global g;
        TF1* mRadiative;
};

extern double gRadiaFcn(double* x, double* p);

extern radiative* gRadiative;

#endif
