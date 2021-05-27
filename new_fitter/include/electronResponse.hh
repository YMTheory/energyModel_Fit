#ifndef _electronResponse_h
#define _electronResponse_h

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
class electronResponse
{
    public:
        electronResponse();
        ~electronResponse();

        static double getElecNonl(double Etrue);
        
        static void loadSimElecNonl();
        static void loadMinSimElecNonl();
        static void loadElecResol();
        static void EmpiricalFit();

        static void Plot();

        static void FitPlot();

        static void FuncConstruct();

    private:
        static const int m_nSimData = 809;
        static double m_SimEtrue[m_nSimData];
        static double m_SimNonl[m_nSimData];

        static double m_scale;

        static bool m_loadSimFile;
        static bool m_doFit;
        static bool m_loadResol;

    private:
        static TGraph* gSimData;
        static TF1* fElecNonl;

    private:
        static double m_p0;
        static double m_p1;
        static double m_p2;
        static double m_p3;

        static double m_q0;
        static double m_q1;
        static double m_q2;

        static double m_c0;
        static double m_c1;
        static double m_c2;
        static double m_s0;

        static double m_ra;
        static double m_rb;
        static double m_rc;

    public:
        static double getp0()       {return m_p0;}
        static void setp0(double p0){m_p0 = p0;}
        static double getp1()       {return m_p1;}
        static void setp1(double p1){m_p1 = p1;}
        static double getp2()       {return m_p2;}
        static void setp2(double p2){m_p2 = p2;}
        static double getp3()       {return m_p3;}
        static void setp3(double p3){m_p3 = p3;}
        static bool getLoadResol()  {return m_loadResol;}

        static double getq0()       {return m_q0;}
        static void setq0(double q0){m_q0 = q0;}
        static double getq1()       {return m_q1;}
        static void setq1(double q1){m_q1 = q1;}
        static double getq2()       {return m_q2;}
        static void setq2(double q2){m_q2 = q2;}

        static double getra()       {return m_ra;}
        static void setra(double ra){m_ra = ra;}
        static double getrb()       {return m_rb;}
        static void setrb(double rb){m_rb = rb;}
        static double getrc()       {return m_rc;}
        static void setrc(double rc){m_rc = rc;}

        static double getc0()       {return m_c0;}
        static void setc0(double c0){m_c0 = c0;}
        static double getc1()       {return m_c1;}
        static void setc1(double c1){m_c1 = c1;}
        static double getc2()       {return m_c2;}
        static void setc2(double c2){m_c2 = c2;}
        static double gets0()       {return m_s0;}
        static void sets0(double s0){m_s0 = s0;}

        static void SetParameters();
        static double calcElecNonl(double E)   {return fElecNonl->Eval(E);}

        static TGraphErrors* gMinElecNonl;
        static TGraphErrors* gElecResol;
        static TF1* fElecResol;
        static TF1* fCerPESigma;
        static TF1* fSctPESigma;
        static TF1* fNPESigma;

};

extern double gElecResolFunc(double* x, double* p);
extern double gCerPESigma(double *x, double *p);
extern double gSctPESigma(double *x, double *p);
extern double gNPESigmaFunc(double* x, double* p);

#endif
