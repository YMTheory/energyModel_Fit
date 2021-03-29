#ifndef _electronResponse_h
#define _electronResponse_h

#include <TGraph.h>
#include <TF1.h>
class electronResponse
{
    public:
        electronResponse();
        ~electronResponse();

        static double getElecNonl(double Etrue);
        
        static void loadSimElecNonl();

        static void Plot();

        static void EmpiricalFit();

    private:
        static const int m_nSimData = 809;
        static double m_SimEtrue[m_nSimData];
        static double m_SimNonl[m_nSimData];

        static double m_scale;

        static bool m_loadSimFile;
        static bool m_doFit;

    private:
        static TGraph* gSimData;
        static TF1* func;

};

#endif
