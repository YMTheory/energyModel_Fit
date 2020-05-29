#include "global.hh"
#include <TF1.h>

#ifndef screening_h
#define screening_h

class screening {
    public:
        screening();
        double screenCorr(double betaE);
        void Plot();

    private:
        global g;
        TF1* mScreen;
};

extern double gScreenFcn(double* x, double* p);

extern screening* gScreen;

#endif
