#ifndef _junoB12Spectrum_h
#define _junoB12Spectrum_h

#include "junoSpectrum.hh"

#include <string>

using namespace std;

class junoB12Spectrum : public junoSpectrum
{
    public:
        junoB12Spectrum(string name);
        ~junoB12Spectrum();

    public:
        void InitTheo();
};

#endif
