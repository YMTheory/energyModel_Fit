#include "global.hh"

#ifndef B12prediction_h
#define B12prediction_h

class B12prediction {

    public:
        double predSpec(double betaE);
        void Plot();
    private:
        global g;

};
#endif
