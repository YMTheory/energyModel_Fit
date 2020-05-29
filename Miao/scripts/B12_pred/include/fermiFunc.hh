#include "global.hh"
#include <TF1.h>

#ifndef fermiFunc_h
#define fermiFunc_h

class fermiFunc{
    public:
        fermiFunc();
        double phasespace(double betaE);
        double fermiCorr(double betaE);
        void Plot();

    private:
        global g;
        TF1* mFermiFcn;

};
extern double gFermiFcn(double* x , double* p);

extern fermiFunc* gFermiFunc;


#endif
