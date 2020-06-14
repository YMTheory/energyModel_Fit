#ifndef _gammaResol_h
#define _gammaResol_h

#include <string>
#include <vector>
#include <TH1.h>

class gammaResol
{
    public:
        gammaResol(std::string name ,
                   double minPE,
                   double maxPE,
                   int nbins
                   );
        ~gammaResol();

        
    public:
        void calcGammaNPE();
        //void predGammaNPE();
        double GetEvis          () {return m_Evis;}
        double GetResolPred     () {return m_resCalc;}
        double GetResolData     () {return m_resData;}
        double GetResolDataErr  () {return m_resDataErr;}
        void Plot();
        double GetChi2();

    private:
        void LoadElecNLData  ();
        void LoadResData     ();

    private:
        std::string m_name;
        double m_minPE;
        double m_maxPE;
        int m_nbins;
        
        std::vector<double> EprmElec;
        double m_mean[1000];
        double m_sigma[1000];
        TH1D* h_totalPE;
        int m_nSamples = 10000;
        //double m_sampleTotPE[10000];

        double m_scale = 1481.06;

        bool m_loadNL;
        double elecNonl[799];
        double m_NLResol = 0.01;

        bool m_loadRes;
        double m_Evis;
        double m_resCalc;
        double m_resData;
        double m_resDataErr;


};
#endif
