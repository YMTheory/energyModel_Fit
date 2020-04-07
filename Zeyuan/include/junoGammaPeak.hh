#ifndef junoGammaPeak_h
#define junoGammaPeak_h

#include <iostream>
#include <string>
#include "TFile.h"
#include "TGraph.h"
#include "junoParameters.hh"
#include "junoEnergyModel.hh"

using namespace std;

class junoGammaPeak
{
    public:
        junoGammaPeak();
        junoGammaPeak(string peakName,
                     string pdfName,
                     double eTru_total,
                     double eTru_single);
        ~junoGammaPeak();

        void   Init(string peakName,  string pdfName,
                    double eTru_total,double eTru_single);
        void   SetEVis       (double val) ;
        void   SetEVisError  (double val) ;
        void   UpdateTheoNL  ();          
        void   UpdateDataNL  ();  

        double GetChi2       ();        

        string GetName       () { return m_name; }
        double GetEVis       () { return m_eVis; }
        double GetEVisError  () { return m_eVisError; }
        double GetETruSingle () { return m_eTru_single; }
        double GetETruTotal  () { return m_eTru_total; }
        double GetTheoScintNL() {return m_theoScintNL;}
        double GetDataScintNL() {return m_dataScintNL;}
        


    private:
        static int s_count;
        string m_name;
        double m_eTru_single; 
        double m_eTru_total; 
        double m_eVis; 
        double m_eVisError; 
        double m_theoScintNL;
        double m_dataScintNL;

        static const unsigned int m_nMaxPdf = 1000;  
        unsigned int m_nPdf;  
        double        m_pdf_eTru [m_nMaxPdf];
        double        m_pdf_prob [m_nMaxPdf];


};

#endif
