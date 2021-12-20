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
        static double getNewAnaCerPE(double E);

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
        static double getA1()            {return m_A1;}
        static double getA2()            {return m_A2;}
        static double getA3()            {return m_A3;}
        static double getA4()            {return m_A4;}
        static void setA1(double A1)     {m_A1 = A1;}
        static void setA2(double A2)     {m_A2 = A2;}
        static void setA3(double A3)     {m_A3 = A3;}
        static void setA4(double A4)     {m_A4 = A4;}

    private:
        static double m_p0;
        static double m_p1;
        static double m_p2;
        static double m_p3;
        static double m_p4;

    public:
        static double getp0()                   {return m_p0;}
        static double getp1()                   {return m_p1;}
        static double getp2()                   {return m_p2;}
        static double getp3()                   {return m_p3;}
        static double getp4()                   {return m_p4;}
        static void setp0(double p0)            {m_p0 = p0;}
        static void setp1(double p1)            {m_p1 = p1;}
        static void setp2(double p2)            {m_p2 = p2;}
        static void setp3(double p3)            {m_p3 = p3;}
        static void setp4(double p4)            {m_p4 = p4;}



};

#endif
