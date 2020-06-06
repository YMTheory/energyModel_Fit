#ifndef junoSpectrum_h
#define junoSpectrum_h

#include <TH1.h>

using namespace std;

class junoSpectrum
{
    public:

        junoSpectrum(int nMaxBins    ,
                    int nMaxBinsData,  
                    int nMaxBr      ,
                    int nMaxGam     ); 
       ~junoSpectrum();

        void    LoadData  (string fileName);
        double  GetChi2   (int nDoF = 0);

        double EvisGamma  (string eTru);
        static double getGammaScale()                    {return m_gammaScale;}
        static void setGammaScale  (double gammaScale)   { m_gammaScale = gammaScale; }


        static double s_n12Ratio;
        TH1F m_eTruH;
        TH1F m_eVisH;
        TH1F m_eSmrH;
        TH1F m_eRecH;

    protected:
        /// junoSpectrum bin number        
        /// number of decay branches		  
        unsigned int m_nBins;       /// 
        unsigned int m_nBinsData;   ///
        unsigned int m_nBr;         /// number of decay branches
        unsigned int m_nGam;        /// max number of gammas per DB

        /// gamma energies      
        double** m_eTruGam;
        double*  m_eTruAlp;
        int** m_eTruGamStr;
        /// continuous beta spectra
        double*  m_binCenter;
        double** m_eTru   ;
        double** m_eTruBck;
        double*  m_eVis   ;
        double*  m_eVisBck;
        double*  m_eRec   ;
        double*  m_eRecBck;
        double*  m_eTheo   ;
        double*  m_eData   ;
        double*  m_eDataErr;


        double m_eMinDraw;
        double m_eMin;
        double m_eMax;
        double m_fitMin;
        double m_fitMax;
        double m_binWidth;
        bool   m_opt;
        bool   m_dataIsLoaded;
        
        int    m_fitMinBin;
        int    m_fitMaxBin;

        string m_name;
        string m_title;


        void DataHist                   (string fileName);
        void TheoHistTree               (string fileName,int isotope=0);
        virtual void InitTheo           ()                 = 0;
        virtual void InitData           (string fileName)  = 0;
        void Normalize                  ();
        void ApplyScintillatorNL        ();
        double GetEVisGamma             (double eTruGam);

        int m_nData;
    
        bool is_positron;
        static double m_pdf_eTru[2000];  // for annihilation e+
        static double m_pdf_prob[2000];
        static double m_max_eTru;

        static double m_gammaScale;


};


#endif
