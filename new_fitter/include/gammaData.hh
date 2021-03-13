#ifndef _gammaData_h
#define _gammaData_h

#include <string>

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
        static double m_scale;

    private:
        static const unsigned int m_nMaxPdf = 2000;
        static double m_pdf_eTrue[m_nMaxPdf];
        static double m_pdf_prob[m_nMaxPdf];
        static double m_max_eTrue;


    private:
        static std::string m_calcOption;  // prmelec ; twolayer

};

#endif
