#ifndef _ELECTRONCERENKOV_WATER_H
#define _ELECTRONCERENKOV_WATER_H

#include <TF1.h>

using namespace std;

class electronCerenkov_water
{
    public:
        electronCerenkov_water();
        ~electronCerenkov_water();

    public:
        // nonlinearity
        static double getCerenkovPE(double E);

        // resolution
        static double getCerenkovPESigma(double N);


    private:
        // nonlinearity
        static double m_w1; 
        static double m_w2; 
        static double m_w3; 
        static double m_w4; 

        // resolution
        static double m_r0; 
        static double m_r1; 
        static double m_r2; 

        static TF1* fCerPE_water;
        static TF1* fCerPESigma_water;



    public:

        static double getw1()             { return m_w1;}
        static double getw2()             { return m_w2;}
        static double getw3()             { return m_w3;}
        static double getw4()             { return m_w4;}
        static void   setw1(double w1)    { m_w1 = w1; }
        static void   setw2(double w2)    { m_w2 = w2; }
        static void   setw3(double w3)    { m_w3 = w3; }
        static void   setw4(double w4)    { m_w4 = w4; }

        static double getr0()             { return m_r0;}
        static double getr1()             { return m_r1;}
        static double getr2()             { return m_r2;}
        static void   setr0(double r0)    { m_r0 = r0; }
        static void   setr1(double r1)    { m_r1 = r1; }
        static void   setr2(double r2)    { m_r2 = r2; }

};


extern double gCerPE_water(double *x, double* p);
extern double gCerPESigma_water(double *x, double* p);


#endif
