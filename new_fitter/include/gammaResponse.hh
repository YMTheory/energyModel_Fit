#ifndef gammaResponse_h
#define gammaResponse_h

#include <string>
#include <TH2D.h>

using namespace std;

class gammaResponse
{
    public:
        gammaResponse(string name, int nBins, double peMin, double peMax);
        ~gammaResponse();

    public:
        void LoadData();
        void LoadPrmElec();
        
        double SampleGamEnergy(int index);
        void calcGamResponse();

        double GetChi2();

    private:
        string m_name;
        int m_nBins;
        double m_peMin;
        double m_peMax;
        double m_Etrue;
        double m_totpeData;
        double m_nonlData;
        double m_nonlDataErr;
        double m_nonlCalc;
        double m_Evis;
        double m_resData;
        double m_resDataErr;

        bool m_loadData;
        bool m_loadPrmElec;

        const int m_nEvents  = 5000;
        const int m_nSamples = 1000;
        TH2D* hPrmElec;

};


#endif
