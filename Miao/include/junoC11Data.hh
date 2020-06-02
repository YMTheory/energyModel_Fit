#ifndef junoC11Data_h
#define junoC11Data_h

#include "junoSpectrum.hh"
#include <string>
using namespace std;

class junoC11Data : public junoSpectrum
{
    public:
        junoC11Data() : junoSpectrum (s_nMaxBins, 
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
