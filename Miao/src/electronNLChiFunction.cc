#include "electronNLChiFunction.hh"
#include "electronNLExperiment.hh"

#include <TMinuit.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>

#include <vector>
#include <iostream>

using namespace std;
double electronNLChiFunction::errElectron = 1.0;  // 100% error 
double electronNLChiFunction::m_chi2 = 0.;
double electronNLChiFunction::m_chi2Min = 10000;

int electronNLChiFunction::m_nParameter = 3;
double electronNLChiFunction::m_bestFit[20] = {0.};
double electronNLChiFunction::m_bestFitError[20] = {0.};

bool electronNLChiFunction::m_DoFit = false;


electronNLChiFunction::electronNLChiFunction()
{;}

electronNLChiFunction::~electronNLChiFunction()
{;}

double electronNLChiFunction::GetChi2(double maxChi2 ) {
    m_chi2 = 0;
    m_chi2 += electronNLExperiment::GetChi2(0);
    if(maxChi2>0 and m_chi2>maxChi2) return maxChi2;

    return m_chi2;
}

void electronNLChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}

double electronNLChiFunction::GetChiSquare( double maxChi2 )
{
    electronNLMinuit = new TMinuit();
    electronNLMinuit->SetFCN(ChisqFCN);
    electronNLMinuit->SetPrintLevel(1);
    
    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    electronNLMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // Configurate parameters
    electronNLMinuit->mnparm(iPar, "kA", 0.96, 0.001, 0.90, 1.00, ierrflag);           iPar++;
    electronNLMinuit->mnparm(iPar, "kB", 5.7e-3, 1e-4, 5.1e-3, 7.0e-3, ierrflag);   iPar++;
    electronNLMinuit->mnparm(iPar, "kC", 1.0, 0.001, 0.90, 1.10, ierrflag);            iPar++;
    //electronNLMinuit->mnparm(iPar, "errElectron", 0.0,  0.01, 0, 0, ierrflag); iPar++;
    
    //electronNLMinuit->FixParameter(0);
    //electronNLMinuit->FixParameter(1);
    //electronNLMinuit->FixParameter(2);

	m_nParameter = electronNLMinuit->GetNumPars();

    // Minimization strategy
    electronNLMinuit->SetErrorDef(1);
    arglist[0]=2;
    electronNLMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 5000; //maxCalls
    arglist[1] = 0.1; // tolerance
    electronNLMinuit->mnexcm("MIGrad", arglist, 2, ierrflag);

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    electronNLMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    // parameters scan ...
    //electronNLMinuit->SetErrorDef(1);
    //TGraph* scan1 = (TGraph*)electronNLMinuit->Contour(100,0,1);
    //TCanvas* cc = new TCanvas();
    //scan1->Draw("APL");
    //cc->SaveAs("tmp.png");

    //arglist[0] = 0;
    //electronNLMinuit->mnexcm("SCAN", arglist, 1, ierrflag);

	for(int i=0; i<m_nParameter; i++)
	{
	    electronNLMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
		//cout<<"curvalue: "<<curvalue<<"	curerror: "<<curerror<<endl;
	}


    m_DoFit = true;

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;
    delete electronNLMinuit;
    m_chi2Min = min;
    return min;
}

void electronNLChiFunction::SetParameters(double *par)
{
    electronQuench::setkA   (par[0]);
    electronQuench::setBirk1(par[1]);
    electronCerenkov::setkC (par[2]);

    //std::cout << par[0] << " " << par[1] << " " << par[2] << endl;
}


void electronNLChiFunction::Plot ()
{
    if(!m_DoFit) {
        cout << " >>> Have Not Done Fitting Yet ! <<< " << endl; return;
    }
    int nPoints = m_nParameter;
    double par[nPoints];
    for(int iPoint=0; iPoint<nPoints; iPoint++) {
        par[iPoint] = m_bestFit[iPoint];
    }
    SetParameters(par);
    electronNLExperiment::Plot();
    
}


void electronNLChiFunction::DrawContour( unsigned int N1, unsigned int N2 )  {
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
            double p0 = m_bestFit[0];  double p1 = 4e-3; double p2 = 1;
            double p1_step = 1e-5;  double p2_step = 0.0001;
            while(p1<8e-3){
                p2 = 1.01;
                while(p2<1.09){
                    electronQuench::setkA(p0);
                    electronQuench::setBirk1(p1);
                    electronCerenkov::setkC(p2);
                    chi2 = GetChi2();
                    p2 += p2_step;
                    cout << p1 << " " << p2 << " " << chi2 << endl;
                } p1 += p1_step;                         
            }
        }

        if(N1 != 1 and N2 != 1) {   // draw kB+kC contour ...
            //TH2F* contour_p1p2 = new TH2F("contour_p1p2", "", );
            double chi2;
            double p0 = 0.970;  double p1 = m_bestFit[1]; double p2 = 1;
            double p0_step = 0.00001;  double p2_step = 0.0001;
            while(p0<0.976){
                p2 = 1.01;
                while(p2<1.09){
                    electronQuench::setkA(p0);
                    electronQuench::setBirk1(p1);
                    electronCerenkov::setkC(p2);
                    chi2 = GetChi2();
                    p2 += p2_step;
                    cout << p0 << " " << p2 << " " << chi2 << endl;
                } p0 += p0_step;                         
            }
        }

        if(N1 != 2 and N2 != 2) {   // draw kA+kB contour ...
            //TH2F* contour_p1p2 = new TH2F("contour_p1p2", "", );
            double chi2;
            double p0 = 0.970;  double p1= 4e-3; double p2 = m_bestFit[2];
            double p0_step = 0.00001;  double p1_step = 0.00001;
            while(p0<0.976){
                p1 = 4e-3;
                while(p1<8e-3){
                    electronQuench::setkA(p0);
                    electronQuench::setBirk1(p1);
                    electronCerenkov::setkC(p2);
                    chi2 = GetChi2();
                    p1 += p1_step;
                    cout << p0 << " " << p1 << " " << chi2 << endl;
                } p0 += p0_step;                         
            }
        }

    }
}
