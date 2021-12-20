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

        static void SaveElecNonl(double* par);

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
        static double m_d0;
        static double m_d1;
        static double m_d2;

        static double m_ra;
        static double m_rb;
        static double m_rc;

        static double m_ma;
        static double m_mb;
        static double m_mc;

        static double m_na;
        static double m_nb;
        static double m_nc;
        static double m_na1;
        static double m_nc1;

        static double m_n1;
        static double m_n2;


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

        static double getma()       {return m_ma;}
        static void setma(double ma){m_ma = ma;}
        static double getmb()       {return m_mb;}
        static void setmb(double mb){m_mb = mb;}
        static double getmc()       {return m_mc;}
        static void setmc(double mc){m_mc = mc;}

        static double getna()       {return m_na;}
        static void setna(double na){m_na = na;}
        static double getnb()       {return m_nb;}
        static void setnb(double nb){m_nb = nb;}
        static double getnc()       {return m_nc;}
        static void setnc(double nc){m_nc = nc;}
        static double getna1()       {return m_na1;}
        static void setna1(double na1){m_na1 = na1;}
        static double getnc1()       {return m_nc1;}
        static void setnc1(double nc1){m_nc1 = nc1;}

        static double getc0()       {return m_c0;}
        static void setc0(double c0){m_c0 = c0;}
        static double getc1()       {return m_c1;}
        static void setc1(double c1){m_c1 = c1;}
        static double getc2()       {return m_c2;}
        static void setc2(double c2){m_c2 = c2;}
        static double gets0()       {return m_s0;}
        static void sets0(double s0){m_s0 = s0;}
        static double getd0()       {return m_d0;}
        static void setd0(double d0){m_d0 = d0;}
        static double getd1()       {return m_d1;}
        static void setd1(double d1){m_d1 = d1;}
        static double getd2()       {return m_d2;}
        static void setd2(double d2){m_d2 = d2;}
        static double getn1()       {return m_n1;}
        static void setn1(double n1){m_n1 = n1;}
        static double getn2()       {return m_n2;}
        static void setn2(double n2){m_n2 = n2;}

        static void SetParameters();
        static double calcElecNonl(double E)   {return fElecNonl->Eval(E);}

        static TGraphErrors* gMinElecNonl;
        static TGraphErrors* gMinElecNPE;
        static TGraphErrors* gElecResol;
        static TF1* fElecResol;
        static TF1* fCerPESigma;
        static TF1* fSctPESigma;
        static TF1* fNPESigma;
        static TF1* fEvisSigma;
        static TF1* fNtotCov;
        static TF1* fEvisCorr;
        static TF1* fEvisNew;
        static TF1* fEvisNew1;

};

extern double gElecResolFunc(double* x, double* p);
extern double gCerPESigma(double *x, double *p);
extern double gSctPESigma(double *x, double *p);
extern double gNtotCov(double* x, double* p);
extern double gNPESigmaFunc(double* x, double* p);
extern double gEvisCorr(double* x, double* p);
extern double gEvisNew(double* x, double* p);
extern double gEvisNew1(double* x, double* p);

#endif
