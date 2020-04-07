#ifndef junoGammaData_h
#define junoGammaData_h

#include "junoData.hh"
#include "junoGammaPeak.hh"
#include "junoParameters.hh"

#include <vector>
#include <iostream>
#include <fstream>
#include "TGraphErrors.h"

class junoGammaData : public junoData
{
    public:
        junoGammaData();

        void LoadData(string fileName);
        double GetChi2( int nDoF = 0 );
        void         GenToyMC();

        void Plot(bool writeToFile);

        TGraphErrors PlotDataScintNL(){return PlotPeaks(1);}
        TGraphErrors PlotTheoScintNL(){return PlotPeaks(2);}


    private:
        vector<junoGammaPeak> m_data;
        TGraphErrors PlotPeaks(int type);
        void AddPeak(string name,string pdfName,
                     double eTruSingle,double eTruTotal);

};

#endif
