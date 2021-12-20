#ifndef michelElectron_water_h
#define michelElectron_water_h

#include <TF1.h>
#include <TGraph.h>

using namespace std;

class michelElectron_water
{

    public:
        michelElectron_water();
        ~michelElectron_water();

    public:
        
        void Initialize();
        void LoadDataSpec();
        void LoadTheoSpec();
        void LoadSimNPE();

        void ApplyLSResponse();

        void Normalize();
        double GetChi2();

        void Plot();

        void SetParameters();

        void Compare();

    public:
        double getc0()           { return m_c0;} 
        void   setc0(double c0)  { m_c0 = c0; }
        double getp1()           { return m_p1;}
        void   setp1(double p1)  { m_p1 = p1;}
        double getp2()           { return m_p2;}
        void   setp2(double p2)  { m_p2 = p2;}


    private:
        bool m_loadData;
        bool m_loadTheo;

        double m_c0;
        double m_p1;
        double m_p2;

        int m_nBin;
        int m_nBinData;
        double m_peMin;
        double m_peMax;
        double m_peMinFit;
        double m_peMaxFit;
        double m_eMin;
        double m_eMax;

    private:
        double m_eBinWidth;
        double m_peBinWidth;
        double m_binCenter[7000];
        double m_eData[100];
        double m_eDataErr[100];
        double m_eTru[7000];
        double m_eVis[7000];
        double m_eTheo[100];


    private:
        TF1* fCerNPE_water;
        TF1* fCerPESigma_water;
        TF1* gaus;

        TGraph* gCerNPE_water;

};

#endif
