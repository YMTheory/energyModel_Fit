#ifndef junoEnergyModel_h
#define junoEnergyModel_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"

using namespace std;

class junoEnergyModel{

    public:
        junoEnergyModel();
        ~junoEnergyModel();

        static double GetKB ()    {return s_kB;}
        static void SetScintP0 (double val) { s_p0 = val; }
        static void SetScintP1 (double val) {s_p1 = val; s_kB = val;
            if(val>99e-4) {s_p1=s_kBMax; s_kB = s_kBMax;}
            if(val<1e-4)  {s_p1=s_kBMin; s_kB=s_kBMin;}
            Update();
        }

		static void SetScintP2 (double val){ s_p2=val;s_cer=val;}

		static void   Load  ();
		static void   Update();
		static double ScintillatorNL(double eTrue);
        static double QuenchNL      (double eTrue);
        static double CerenkovNL    (double eTrue);

		static void SaveCurves();
		static TGraph DrawElectronScintNL(int nSamples = 798,double eMax=7.98);
		static TGraph DrawElectronQuenchNL(int nSamples = 798,double eMax=7.98);
		static TGraph DrawElectronCerenkovNL(int nSamples = 798,double eMax=7.98);
		static vector<double> SampleElectronScintNL  (int nSamples,double eMax=8);
		static vector<double> SampleElectronQuenchNL  (int nSamples,double eMax=8);
		static vector<double> SampleElectronCerenkovNL  (int nSamples,double eMax=8);

		static double  s_cer;
		static double  s_rad;

		static double  s_p0;
		static double  s_p1;
		static double  s_p2;
		static double  s_p3;

	private: 
		static double  s_kB;
		static double  s_kBResid;
		static bool    s_isLoaded;
		static constexpr double s_kBMax = 99e-4;
		static constexpr double s_kBMin = 1e-4;
		static const unsigned int s_nKb      =  100;
		static const unsigned int s_nSamples = 1000;  // 799 energy bins
		static constexpr double s_samplingResol  = 0.01;

		static double s_energySamples          [s_nSamples];
		static double s_quenchingShape1 [s_nKb][s_nSamples];
		static double s_cerenkovShape          [s_nSamples];

		static double* s_quenchingShape1_lowKb;
		static double* s_quenchingShape1_higKb;

		static double PhysicsScintillator    (double eTru);

		static double ScintillatorShape      (double eTru);

		static vector<double> m_energySamples;
};
#endif
