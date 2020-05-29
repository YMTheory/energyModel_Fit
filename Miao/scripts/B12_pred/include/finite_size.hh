#include <global.hh>

#include <TF1.h>

#ifndef finite_size_h
#define finite_size_h

class finite_size{

    public:
        finite_size();
        double elecmagCorr(double betaE);
        double weakCorr(double betaE);
        void Plot();

    private:
        global g;
        TF1* mFiniteSize;
};

extern double gFiniteSizeFcn(double* x, double* p);

extern finite_size* gFiniteSize;

#endif
