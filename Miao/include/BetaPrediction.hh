#ifndef BetaPrediction_h
#define BetaPrediction_h

#include <TMath.h>
#include <TF1.h>

class BetaPrediction {
    public:
        BetaPrediction();

        static double predB12Spec(double betaE);
        static void Plot();
        static void setK(double m_K) { K = m_K; }
        static double getK() {return K;}
        static TF1* B12Fcn;


    private:
        static constexpr double alpha = 1/137.036;
        static constexpr double Z = 12; 
        static constexpr double A = 5;
        static double gamma;
        static double R ;
        static constexpr double E0 = 13.4; //MeV
        static constexpr double W0 = E0/0.511+1;  //MeV
        static double K;

        static const int m_bin = 200;   // 200 energy bins

};

extern double gBetaSpecFcn (double* x, double* p) ;

extern BetaPrediction* gBetaPredicion;

#endif
