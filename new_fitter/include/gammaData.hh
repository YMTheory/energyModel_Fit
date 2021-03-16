#ifndef _gammaData_h
#define _gammaData_h

#include <string>
#include <vector>
#include <TH2D.h>

class gammaData
{
    public:
        gammaData(std::string name,
                  double minPE,
                  double maxPE,
                  int nbins
                 );

        ~gammaData();

    
    public:
        double GetEtrue         () {return m_Etrue;}
        double GetEvis          () {return m_Evis;}
        double GetNonlPred      () {return m_nonlCalc;}
        double GetNonlData      () {return m_nonlData;}
        double GetNonlDataErr   () {return m_nonlDataErr;}
        double GetResolPred     () {return m_resCalc;}
        double GetResolData     () {return m_resData;}
        double GetResolDataErr  () {return m_resDataErr;}
        void LoadData           ();
        void Plot               ();
        double GetChi2          ();

        void LoadGammaData ();
        void LoadPrimElecDist();
        void LoadPrimElecSamples();

        void calcGammaResponse();

    
    private:

        std::string m_name;
        double m_minPE;
        double m_maxPE;
        int m_nbins;

        bool m_loadData;

        double m_Etrue;
        double m_Evis;
        double m_nonlCalc;
        double m_nonlData;
        double m_nonlDataErr;
        double m_resCalc;
        double m_resData;
        double m_resDataErr;

    private:
        double m_scale;

    private:
        static const unsigned int m_nMaxPdf = 600;
        double m_pdf_eTrue[m_nMaxPdf];
        double m_pdf_prob[m_nMaxPdf];
        double m_max_eTrue;

    private:
        TH2D* elec_hist;
        static const unsigned int m_nSamples = 5000;
        double m_mean[m_nSamples];

    private:
        static std::string m_calcOption;  // prmelec ; twolayer

};

#endif
