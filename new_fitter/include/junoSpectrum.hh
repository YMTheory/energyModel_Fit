#ifndef _junoSpectrum_h
#define _junoSpectrum_h

#include <iostream>
#include <string>

using namespace std;

class junoSpectrum
{

    public:
        junoSpectrum()  {;}
        ~junoSpectrum() {;}


    public:
        void InitTheo();
        void TheoHistTree(string theofile);
    //     void InitData();

    
    public:
        string m_name;
        bool m_isPositron;

        int m_nBranch;
        int m_nGamma[10];
        double m_eTruGam[10][10];
        string m_eTruGamStr[10][10];
};

#endif
