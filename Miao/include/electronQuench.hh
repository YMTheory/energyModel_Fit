#ifndef _ELECTRONQUENCH_H
#define _ELECTRONQUENCH_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include "junoParameters.hh"

using namespace std;

class electronQuench
{
    public:
        electronQuench();
        ~electronQuench();
    
    public:
        static void setBirk1(double birk1) { m_birk1 = birk1;}
        static double getBirk1()           {return m_birk1;}
        static void setkA(double val)      { m_kA = val; }
        static double getkA()              {return m_kA;}

        static void LoadStopPowData();
        static double Integral_BirkLaw(double E);

        static void LoadNLData();
        static void Update();

        static double ScintillatorNL (double eTrue);

        static void Plot ();


    private:
        static double m_birk1;
        static double m_kA;
		static const unsigned int m_nKb          = 100;
		static const unsigned int m_nSamples     = 1000;
		static constexpr double m_samplingResol  = 0.05;   // bining: 20keV/bin
		static double  m_kBResid;
        static double p0 ;
        static double p1 ;
        static double p2 ;
        static double p3 ;

        static double m_edep[1000];
        static double m_nonl[1000];

		static double m_quenchingShape1 [m_nKb][m_nSamples];
		static double* m_quenchingShape1_lowKb;
		static double* m_quenchingShape1_higKb;

        static bool m_loadStopPowData;
        static bool m_loadNLData;

        static std::vector<double> m_Etrue;
        static std::vector<double> m_StopPow;

        static double ScintillatorShape    (double eTrue);
        static double SimulationNLShape    (double eTrue);
        static double IntegralNLShape      (double eTrue);
        static double EmpiricalNLShape     (double eTrue);

};

#endif
