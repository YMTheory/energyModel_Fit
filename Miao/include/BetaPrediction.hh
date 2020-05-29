#ifndef BetaPrediction_h
#define BetaPrediction_h

#include <TMath.h>

class BetaPrediction {
    public:
        static double predSpec(double betaE);
        static void Plot();

    private:
        static constexpr double alpha = 1/137.036;
        static constexpr double Z = 12; 
        static constexpr double A = 5;
        static double gamma;
        static double R ;
        static constexpr double E0 = 13.4; //MeV
        static constexpr double W0 = E0/0.511+1;  //MeV
        static constexpr double K = 1/84355.9;

        static const int m_bin = 200;   // 200 energy bins

};

#endif
