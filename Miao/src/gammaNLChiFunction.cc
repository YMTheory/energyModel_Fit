#include "gammaNLChiFunction.hh"
#include "gammaNLExperiment.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

#include <iostream>

using namespace std;

double gammaNLChiFunction::errGamma = 1.0;  // 100% error 
double gammaNLChiFunction::m_chi2 = 0.;

int gammaNLChiFunction::m_nParameter = 3;
double gammaNLChiFunction::m_bestFit[20] = {0.};
double gammaNLChiFunction::m_bestFitError[20] = {0.};

bool gammaNLChiFunction::m_DoFit = false;

gammaNLChiFunction::gammaNLChiFunction()
{;}

gammaNLChiFunction::~gammaNLChiFunction()
{;}

double gammaNLChiFunction::GetChi2 ( double maxChi2 ) {
    m_chi2 = 0;
    m_chi2 += gammaNLExperiment::GetChi2(0);
    m_chi2 += (gammaNLExperiment::m_gammaScale/junoParameters::m_gammaError, 2);
    if(maxChi2>0 and m_chi2>maxChi2) return maxChi2;

    return m_chi2;
}


void gammaNLChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


double gammaNLChiFunction::GetChiSquare(double maxChi2)
{
    gammaNLMinuit = new TMinuit();
    gammaNLMinuit->SetFCN(ChisqFCN);
    gammaNLMinuit->SetPrintLevel(1);
    
    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    gammaNLMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // Configurate parameters
    gammaNLMinuit->mnparm(iPar, "kA", 0.96, 0.01, 0.9, 1.0, ierrflag); iPar++;
    gammaNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-4, 6e-3, 7e-3, ierrflag); iPar++;
    gammaNLMinuit->mnparm(iPar, "kC", 1.0, 0.01, 0.9, 1.1, ierrflag); iPar++;
    gammaNLMinuit->mnparm(iPar, "errGamma", 0, 0.01, 0, 0, ierrflag); iPar++;
    
    //gammaNLMinuit->FixParameter(0);
    //gammaNLMinuit->FixParameter(1);
    //gammaNLMinuit->FixParameter(2);
    //gammaNLMinuit->FixParameter(3);

    // Minimization strategy
    gammaNLMinuit->SetErrorDef(1);
    arglist[0]=2;
    gammaNLMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 5000; //maxCalls
    arglist[1] = 0.01; // tolerance
    gammaNLMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    gammaNLMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

	for(int i=0; i<m_nParameter; i++)
	{
	    gammaNLMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
		//cout<<"curvalue: "<<curvalue<<"	curerror: "<<curerror<<endl;
	}

    m_DoFit = true;

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;
    delete gammaNLMinuit;
    return min;
}


void gammaNLChiFunction::SetParameters(double *par) {
    // 4-parameter fitting
    electronQuench::setkA             (par[0]);
    electronQuench::setBirk1          (par[1]);
    electronCerenkov::setkC           (par[2]);
    gammaNLExperiment::setGammaScale  (par[3]);
    
    // 3-parameter fitting
    //electronQuench::setBirk1          (par[0]);
    //electronCerenkov::setkC           (par[1]);
    //gammaNLExperiment::setGammaScale  (par[2]);
    //electronQuench::setkA             ((1-par[1]*58.517/1481.06)/0.9796);
}


void gammaNLChiFunction::Plot ()
{
    if(!m_DoFit) {
        cout << " >>> Have Not Done Fitting Yet ! <<< " << endl; return;
    }
    int nPoints = m_nParameter;
    double par[nPoints];
    for(int iPoint=0; iPoint<nPoints; iPoint++) {
        par[iPoint] = m_bestFit[iPoint];
    }
    //double par[3] = {0.967, 6.5e-3, 1};
    //SetParameters(par);
    gammaNLExperiment::Plot();
}

void gammaNLChiFunction::DrawContour( unsigned int N1, unsigned int N2 )  {
    if(!m_DoFit) {
        cout << " >>> Not Fitting Yet ! <<< "; return; 
    } else if(N1==N2) {
        cout << " >>> Error: Two Same Parameter No <<< " << endl; return;
    } else if(N1 >= m_nParameter or N2 >=m_nParameter) {
        cout << " >>> Error: Parameter No Out of Range <<< "; return;
    } else {
        cout << " >>> Draw Contour Plot <<< " << endl;
        if(N1 != 0 and N2 != 0) {   // draw kB+kC contour ...
            //TH2F* contour_p1p2 = new TH2F("contour_p1p2", "", );
            double chi2;
            double p0 = m_bestFit[0];  double p3 = m_bestFit[3];
            double p1 = m_bestFit[1]-3*m_bestFitError[1]; double p2 = m_bestFit[2]-3*m_bestFitError[2];
            double p1_step = 0.00001;  double p2_step = 0.001;
            while(p1<m_bestFit[1]+3*m_bestFitError[1]){
                p2 = m_bestFit[2]-3*m_bestFitError[2];
                while(p2<m_bestFit[2]+3*m_bestFitError[2]){
                    electronQuench::setkA(p0);
                    electronQuench::setBirk1(p1);
                    electronCerenkov::setkC(p2);
                    gammaNLExperiment::setGammaScale  (p3);
                    chi2 = GetChi2();
                    p2 += p2_step;
                    cout << p1 << " " << p2 << " " << chi2 << endl;
                } p1 += p1_step;                         
            }
        }

        if(N1 != 1 and N2 != 1) {   // draw kB+kC contour ...
            //TH2F* contour_p1p2 = new TH2F("contour_p1p2", "", );
            double chi2;
            double p0 = m_bestFit[0]-3*m_bestFitError[0];  double p1 = m_bestFit[1]; 
            double p2 = m_bestFit[2]-3*m_bestFitError[2];  double p3 = m_bestFit[3];
            double p0_step = 0.0001;  double p2_step = 0.001;
            while(p0<m_bestFit[0]+3*m_bestFitError[0]){
                p2 = m_bestFit[2]-3*m_bestFitError[2];
                while(p2<m_bestFit[2]+3*m_bestFitError[2]){
                    electronQuench::setkA(p0);
                    electronQuench::setBirk1(p1);
                    electronCerenkov::setkC(p2);
                    gammaNLExperiment::setGammaScale  (p3);
                    chi2 = GetChi2();
                    p2 += p2_step;
                    cout << p0 << " " << p2 << " " << chi2 << endl;
                } p0 += p0_step;                         
            }
        }

     /*   if(N1 != 2 and N2 != 2) {   // draw kA+kB contour ...
            //TH2F* contour_p1p2 = new TH2F("contour_p1p2", "", );
            double chi2;
            double p0 = m_bestFit[0]-5*m_bestFitError[0];  double p1= m_bestFit[1]-m_bestFitError[1]*5; 
            double p2 = m_bestFit[2]; double p3 = m_bestFit[3];
            double p0_step = 0.0001;  double p1_step = 0.00001;
            while(p0<m_bestFit[0]+5*m_bestFitError[0]){
                p1 = 4e-3;
                while(p1<m_bestFit[1]+5*m_bestFitError[1]){
                    electronQuench::setkA(p0);
                    electronQuench::setBirk1(p1);
                    electronCerenkov::setkC(p2);
                    gammaNLExperiment::setGammaScale  (p3);
                    chi2 = GetChi2();
                    p1 += p1_step;
                    cout << p0 << " " << p1 << " " << chi2 << endl;
                } p0 += p0_step;                         
            }
        } */

        // for two parameters fitting
        if(N1!=2 and N2!=2) {
            double chi2;
            double p0 = m_bestFit[0]-1*m_bestFitError[0];  double p1; double p2 = m_bestFit[2];
            double p0_step = 0.00001; double p1_step = 0.001;
            while(p0<m_bestFit[0]+3*m_bestFitError[0]) {
                p1 = m_bestFit[1]-3*m_bestFitError[1];
                while(p1<m_bestFit[1]+3*m_bestFitError[1]) {
                    electronQuench::setBirk1(p0);
                    electronCerenkov::setkC(p1);
                    electronQuench::setkA(1-p1*58.517/1481.06);
                    gammaNLExperiment::setGammaScale(p2);
                    chi2 = GetChi2();
                    p1 += p1_step;
                    cout << p0 << " " << p1 << " " << chi2 << endl;
                } p0 += p0_step;
            }
        }


    }
}
