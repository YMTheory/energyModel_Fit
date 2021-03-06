#ifndef _ELECTRONCERENKOV_H
#define _ELECTRONCERENKOV_H

#include <vector>
#include "junoParameters.hh"

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

    private:
        static double m_kC;
        static double m_energyScale;

        static bool m_LoadCerenkov;

        static std::vector<double> m_Etrue;
        static std::vector<double> m_Cerenkov;

};

#endif
