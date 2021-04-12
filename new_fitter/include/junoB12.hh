#ifndef _junoB12_h
#define _junoB12_h

#include "gammaData.hh"
#include <TH1.h>

class junoB12
{

    public:
        junoB12();
        ~junoB12();

    public:
        void Initialize();

        void LoadDataSpec();
        void LoadTheoSpec();

        void ApplyScintillatorNL();
        void Normalize();
        
        double GetChi2();

        void Plot();


    private:
        bool m_loadData;
        bool m_loadTheo;

        double m_binWidth;
        double m_binCenter[14000];
        double m_eData[14000];
        double m_eDataErr[14000];
        double m_eTheo[100];
        double m_eTru[3][100];
        double m_eVis[100];
        double m_eTruGam[3][2];
        double m_eVisGam[3][2];

        gammaData* gamma4440;
        gammaData* gamma3215;
};
#endif
