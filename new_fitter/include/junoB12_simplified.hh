#ifndef junoB12_simplified_h
#define junoB12_simplified_h

#include <TF1.h>

using namespace std;

class junoB12_simplified
{
    
    public:
        junoB12_simplified(int nBinsData, double fitMinE, double fitMaxE);
        ~junoB12_simplified();

    public:
        
        void Initialize();
        void LoadDataSpec();
        void LoadTheoSpec();

        void ApplyResponse();

        void Normalize();
        double GetChi2();

        void Plot();

    private:
        bool m_loadData;
        bool m_loadTheo;

        int m_nBin;
        int m_nBinData;
        double m_fitMinPE;
        double m_fitMaxPE;
        double m_eMin;
        double m_eMax;
        double m_peMin;
        double m_peMax;


        double m_eBinWidth;
        double m_peBinWidth;
        double m_eBinCenter[1500];
        double m_peBinCenter[1500];
        double m_peData[200];
        double m_peDataErr[200];
        double m_eTru[1500];
        double m_eVis[1500];
        double m_peTheo[200];

    private:
        TF1* gaus;

};



#endif
