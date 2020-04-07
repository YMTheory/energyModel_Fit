#ifndef _ELECTRONQUENCH_H
#define _ELECTRONQUENCH_H

#include <vector>

using namespace std;

class electronQuench
{
    public:
        electronQuench();
        ~electronQuench();
    
    public:
        void setBirk1(double birk1) { m_birk1 = birk1; }
        double getBirk1() {return m_birk1;}

        double Integral_BirkLaw(double kB1, double E);


    private:
        static double m_birk1;
		static const unsigned int s_nKb      = 256;
		static const unsigned int s_nSamples = 1000;

		static double s_quenchingShape1 [s_nKb][s_nSamples];

        std::vector<double> Etrue;
        std::vector<double> StopPow;

};

#endif
