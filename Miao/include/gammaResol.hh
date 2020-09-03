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
        double GetEtrue         () {return m_Etrue;}
        double GetEvis          () {return m_Evis;}
        double GetNonlPred      () {return m_nonlCalc;}
        double GetNonlData      () {return m_nonlData;}
        double GetNonlDataErr   () {return m_nonlDataErr;}
        double GetResolPred     () {return m_resCalc;}
        double GetResolData     () {return m_resData;}
        double GetResolDataErr  () {return m_resDataErr;}
        void LoadData        ();
        void Plot();
        double GetChi2();
        void check_nonl();

    private:
        void LoadElecNLData  ();
        void LoadGammaNLData ();
        void LoadResData     ();
        void LoadPrimElecDist();

        double interpolate_nonl(int idx, double E);

    private:
        std::string m_name;
        double m_minPE;
        double m_maxPE;
        int m_nbins;
        
        std::vector<double> EprmElec;
        double m_mean[5000];
        double m_sigma[5000];
        TH1D* h_totalPE;
        int m_nSamples = 5000;
        //double m_sampleTotPE[10000];

        double m_scale = 3350/2.22;   // prior energy scale, nH energy scale

        bool m_loadData;
        bool m_FitNL;
        bool m_FitRes;
        unsigned int m_NL_option;

        double elecEtrue[810];
        double elecNonl[810];
        double m_NLResol = 0.01;

        double m_Etrue;
        double m_Evis;
        double m_nonlCalc;
        double m_nonlData;
        double m_nonlDataErr;
        double m_resCalc;
        double m_resData;
        double m_resDataErr;


        // data for primary electron distribution: 
        static const unsigned int m_nMaxPdf = 2000;
        int m_max_eTrue;
        double m_pdf_eTrue[m_nMaxPdf];
        double m_pdf_prob[m_nMaxPdf];



};
#endif
