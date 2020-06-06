#ifndef junoC10Data_h
#define junoC10Data_h

#include "junoSpectrum.hh"
#include <string>
using namespace std;

class junoC10Data : public junoSpectrum
{
    public:
        junoC10Data() : junoSpectrum (s_nMaxBins, 
                                      s_nMaxBinsData, 
                                      s_nMaxBr,
                                      s_nMaxGam) {;}
        void SetParameters           ();
        void InitTheo                ();
        void InitData                (string fileName);
        void Plot                    ();


    private:
        static unsigned int s_nMaxBr;
        static unsigned int s_nMaxGam;
        static unsigned int s_nMaxBins;
        static unsigned int s_nMaxBinsData;

};



#endif
