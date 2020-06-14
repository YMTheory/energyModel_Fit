#ifndef _ELECTRONRESOL_H
#define _ELECTRONRESOL_H

class electronResol
{
    public:
        electronResol();
        ~electronResol();

    public:
        static void setpA(double val) { m_pA = val; }
        static void setpB(double val) { m_pB = val; }
        static void setpC(double val) { m_pC = val; }

        static double getpA()         {return m_pA;}
        static double getpB()         {return m_pB;}
        static double getpC()         {return m_pC;}

        static double energySmear(double Evis);
        static double Resolution(double Evis);

    private:
        static double m_pA;
        static double m_pB;
        static double m_pC;

};
#endif
