#ifndef junoExperiment_h
#define junoExperiment_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

class junoExperiment
{
    public:
        virtual void   LoadData() = 0;
        virtual double GetChi2( int nDoF = 0 ) = 0;
        
    protected:
        int m_nData;
};
#endif
