#ifndef junoB12Data_h
#define junoB12Data_h

#include "junoSpectrum.hh"
#include <string>

using namespace std;

class junoB12Data : public junoSpectrum
{
  public:
    junoB12Data() : junoSpectrum(s_nMaxBins,
                               s_nMaxBinsData,
                               s_nMaxBr,
                               s_nMaxGam){;}
    void  SetParameters();
    void  InitTheo     ();
    void  InitData     (string fileName);
    void setNuWM      (double nuWM) {m_nuWM = nuWM;}
    double getNuWM    () {return m_nuWM;}

    void Plot          ();
    
  private:
    static unsigned int s_nMaxBr;
    static unsigned int s_nMaxGam;
    static unsigned int s_nMaxBins;
    static unsigned int s_nMaxBinsData;

    static double m_nuWM;

};


#endif
