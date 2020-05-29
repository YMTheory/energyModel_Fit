#include <TMath.h>

#ifndef global_h
#define global_h

class global{
    // global variables 
    public:
        double E0 = 13.4;  // end point -> Q value
        double W0 = E0/0.511+1;
        double alpha = 1/137.036;
        double Z = 12;
        double A = 5;
        double gamma = TMath::Sqrt(1-(alpha*Z));
        double R = 0.0029*TMath::Power(A,1/3)+0.0063*TMath::Power(A,-1/3)-0.017*TMath::Power(A,-1);
};

#endif
