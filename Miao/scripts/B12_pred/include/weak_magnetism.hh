#include "global.hh"
#include <TF1.h>

#ifndef weak_magnetism_h
#define weak_magnetism_h

class weak_magnetism{
    public:
        weak_magnetism();
        double weakMagCorr(double betaE);
        void Plot();

    private:
        global g;
        TF1* mWeakMag;
};

#endif
