#ifndef michelElectron_h
#define michelElectron_h

#include <TH1.h>
#include <TF1.h>

using namespace std;

class michelElectron
{
    public:
        michelElectron();
        ~michelElectron();

    public:
        void Initialize();

        void LoadMCTruth();
        void LoadMCData();

        void CalcVisibleEnergy();

        
        double GetChi2();

        void Plot();

    private:

        bool m_loadDeposit;
        bool m_loadVisible;

        int m_nBin;
        int m_nBinData;
        double m_binWidth;
        double m_binDataWidth;
        double m_binCalcWidth;

        double m_energyScale;

        double m_dataEntries;
        double m_calcEntries;

        double m_EdepMin;
        double m_EdepMax;
        double m_EvisMin;
        double m_EvisMax;

    private:
        TH1D* hEdep;
        TH1D* hEvisData;
        TH1D* hEvisCalc;
        TH1D* hTmpEvisCalc;

        TF1* gaus;
};



#endif
