#ifndef _junoSpectrum_h
#define _junoSpectrum_h

#include <iostream>
#include <string>
#include <map>

using namespace std;

class junoSpectrum
{

    public:
        junoSpectrum(int nMaxBins    ,
                     int nMaxBinsData,  
                     int nMaxBr      ,
                     int nMaxGam     ,
                     double eMin     ,
                     double eMax     ,
                     string name    ); 

        ~junoSpectrum() ;


    public:
        void LoadData();
        void InitTheo();
        void TheoHistTree(string theofile);
        void InitData();
        void DataHistTree(string datafile);

        void ApplyScintillatorNL();
        void LoadPrmElecDist();
        double EvisGamma(int Etrue);

        double GetChi2();

    
    public:
        string m_name;
        bool m_isPositron;

    public:
        // decay branch
        int m_nBranch;        // branch number
        int m_nGam;           // Max gamma number among all branches

        // gamma number & energy 
        int* m_nGamma;
        double** m_eTruGam;
        int** m_eTruGamStr;
    
        // beta spectrum
        int m_nBins;
        double m_binWidth;
        double m_eMin;
        double m_eMax;
    
        double m_fitMin;
        double m_fitMax;
        int m_fitMinBin;
        int m_fitMaxBin; 

        int m_nBinsData;
        double** m_eTru;
        double* m_binCenter;
        double* m_eData;
        double* m_eDataErr;

        // alpha energy
        double* m_eTruAlp;

    public:  // nonlinearity apply
        double* m_eVis;
        int m_max_eTru;

        int m_nPdfBins;
        map<int, int> mapPdfMaxEtrue;
        map<int, double*> mapPdfEtrue;
        map<int, double*> mapPdfProb;

};

#endif
