#ifndef _ELECTRONCERENKOV_H
#define _ELECTRONCERENKOV_H

#include <vector>
#include <TGraph.h>

using namespace std;

class electronCerenkov
{
    public:
        electronCerenkov();
        ~electronCerenkov();
    
    public:
        static void setkC(double val)           {m_kC = val;}
        static double getkC()                   {return m_kC;}
        static void setEnergyScale(double val)  {m_energyScale = val;}
        static double getEnergyScale()          {return m_energyScale;}

        static void LoadCerenkov();
        
        static double getCerenkovPE(double E);

        static double getCerPE(double E);

        static double getAnaCerPE(double E);
        static double getSimCerPE(double E);

        static void Plot();

    private:
        static double m_kC;
        static double m_energyScale;

        static bool m_LoadCerenkov;

        static std::vector<double> m_Etrue;
        static std::vector<double> m_Cerenkov;

        static const int m_nData = 900;
        static double m_E[m_nData];
        static double m_nonl[m_nData];

        static TGraph* gNPE_elec;

        static double m_A1;
        static double m_A2;
        static double m_A3;
        static double m_A4;
        static double m_E0;

    public:
        double getA1()            {return m_A1;}
        double getA2()            {return m_A2;}
        double getA3()            {return m_A3;}
        double getA4()            {return m_A4;}
        void setA1(double A1)     {m_A1 = A1;}
        void setA2(double A2)     {m_A2 = A2;}
        void setA3(double A3)     {m_A3 = A3;}
        void setA4(double A4)     {m_A4 = A4;}
};

#endif
