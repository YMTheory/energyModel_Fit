#ifndef electronResExperiment_h
#define electronResExperiment_h

#include <vector>

#include "TH1D.h"
#include "TFile.h"

class electronResExperiment
{
    public:
        electronResExperiment();
        ~electronResExperiment();

    public:
        static void LoadData();
        static void LoadB12Data();
        static void LoadC11Data();

        static void UpdateTrueSpectrum();
        static void UpdateTheoSpectrum();

        static double GetChi2();
        static void   PlotB12();
        static void   PlotC11();

    private:
        static bool m_LoadData;
        static bool m_LoadB12Data;
        static bool m_LoadC11Data;
        static bool m_CalcTheo;

        static double m_B12Energy[500];
        static double m_B12Origin[500];
        static double m_B12Nonl[500];
        static double m_B12Smear[500];
        static double m_B12Fit[500];
        static std::vector<double> B12energySpec;
        static std::vector<double> C11energySpec;

        static TH1D* h_B12Origin;
        static TH1D* h_B12Nonl;
        static TH1D* h_B12Smear;
        static TH1D* h_B12Fit;
        static TFile* B12file;

        static TH1D* h_C11Origin;
        static TH1D* h_C11Nonl;
        static TH1D* h_C11Smear;
        static TH1D* h_C11Fit;
        static TFile* C11file;
};
#endif
