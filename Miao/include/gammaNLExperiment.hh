#ifndef _GAMMANLEXPERIMENT_H
#define _GAMMANLEXPERIMENT_H

#include <vector>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "junoParameters.hh"

using namespace std;

class electronQuench;
class electronCerenkov;
class primaryElectronDistribution;
class gammaNLExperiment
{
    public:
        //gammaNLExperiment(electronQuench* aQuench, electronCerenkov* aCerenkov, primaryElectronDistribution* aPED);
        gammaNLExperiment();
        ~gammaNLExperiment();

    public:

        static double m_gammaScale;
        static void   setGammaScale   (double val)  {m_gammaScale = val;}
        static double getGammaScale ()            {return m_gammaScale;}
        
        static void   LoadData();
        static double GetChi2( int nDoF = 0 );

        static void   LoadPrimaryElecDist();

        void Calculate(double *apar);

        static void   UpdateDataGammaNL();
        static void   UpdateTheoGammaNL();

        static double CalculateGammaNL( int iSource );
        
        static void   Plot       ();


    private:
        static TGraphErrors* mTrueGammaNL;
        static TGraph* mFitGammaNL;

        static bool m_LoadGammaData;
        static bool m_LoadPrmElecDist;
        static bool m_CalcTheo;


        static vector<double> Etrue;
        static vector<string> source_name;

        static const unsigned int m_nMaxPdf = 1000;
        static const unsigned int m_nMaxSources = 20;
        static double m_pdf_eTru[m_nMaxSources][m_nMaxPdf];
        static double m_pdf_prob[m_nMaxSources][m_nMaxPdf];
        static double m_max_eTru[m_nMaxSources];
        static double m_ratio[m_nMaxSources];

        static double m_scale;
};

#endif
