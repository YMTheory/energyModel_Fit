#include "junoNLChiFunction.hh"
#include "electronNLExperiment.hh"
#include "gammaNLExperiment.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "BetaPrediction.hh"
#include "junoB12Data.hh"
#include "junoC11Data.hh"
#include "junoParameters.hh"

#include <iostream>
#include <fstream>

using namespace std;

junoB12Data* junoNLChiFunction::m_b12Data = 0;
junoC11Data* junoNLChiFunction::m_c11Data = 0;
junoC10Data* junoNLChiFunction::m_c10Data = 0;

double junoNLChiFunction::m_chi2    = 0.;
double junoNLChiFunction::m_chi2Min = 100000;
double junoNLChiFunction::m_chi2B12 = 0;
double junoNLChiFunction::m_chi2C11 = 0;
double junoNLChiFunction::m_chi2C10 = 0;
double junoNLChiFunction::m_chi2Gam = 0;

int junoNLChiFunction::m_nParameter = 3;
double junoNLChiFunction::m_bestFit[20] = {0.};
double junoNLChiFunction::m_bestFitError[20] = {0.};

bool junoNLChiFunction::m_DoFit = false;
bool junoNLChiFunction::m_gridSearch = false;

double junoNLChiFunction::final_kA = 0;
double junoNLChiFunction::final_kB = 0;
double junoNLChiFunction::final_kC = 0;
double junoNLChiFunction::m_kALimit = 0.01;
double junoNLChiFunction::m_kBLimit = 3e-4;
double junoNLChiFunction::m_kCLimit = 0.01;

junoNLChiFunction::junoNLChiFunction() {
    m_b12Data = new junoB12Data();
    m_b12Data->SetParameters();
    m_c11Data = new junoC11Data();
    m_c11Data->SetParameters();
    m_c10Data = new junoC10Data();
    m_c10Data->SetParameters();
}

junoNLChiFunction::~junoNLChiFunction() {
    delete m_b12Data;
    delete m_c11Data;
    delete m_c10Data;
}

void junoNLChiFunction::LoadData()
{
	std::cout << " ---> Start loading data " << std::endl;
    m_b12Data   ->LoadData(junoParameters::B12DataFile);
    m_c11Data   ->LoadData(junoParameters::C11DataFile);
    m_c10Data   ->LoadData(junoParameters::C10DataFile);
}


double junoNLChiFunction::GetChi2(double maxChi2) {
    m_chi2    = 0;
    m_chi2Gam = 0;
    m_chi2B12 = 0;
    m_chi2C11 = 0;
    m_chi2C10 = 0;
    //m_chi2 += electronNLExperiment::GetChi2(0);

    if(junoParameters::fitGammaSources) {
        m_chi2Gam = gammaNLExperiment::GetChi2();
        m_chi2 += m_chi2Gam;
    }

    if(junoParameters::fitB12) {
        m_chi2B12 = m_b12Data ->GetChi2();
        m_chi2 += m_chi2B12;
    }

    if(junoParameters::fitC11) {
        m_chi2C11 = m_c11Data ->GetChi2();
        m_chi2 += m_chi2C11;
    }

    if(junoParameters::fitC10) {
        m_chi2C10 = m_c10Data ->GetChi2();
        m_chi2 += m_chi2C10;
    }

    if(maxChi2>0 and m_chi2>maxChi2) return maxChi2;
    return m_chi2;
}


void junoNLChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


void junoNLChiFunction::SetParameters(double* par)
{
    //BetaPrediction::setK                (par[0]);   // B12 Spec Normalization
    electronQuench::setkA               (par[0]);
    electronQuench::setBirk1            (par[1]);
    electronCerenkov::setkC             (par[2]);
    gammaNLExperiment::setGammaScale    (par[3]);
    junoSpectrum::setGammaScale         (par[3]);
    //electronQuench::setkA               ((1-par[2]*58.517/1481.06)/0.9796);
}



double junoNLChiFunction::GetChiSquare(double maxChi2)
{
    junoNLMinuit = new TMinuit();
    junoNLMinuit->SetFCN(ChisqFCN);
    junoNLMinuit->SetPrintLevel(1);
    
    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    junoNLMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // Configurate parameters
    junoNLMinuit->mnparm(iPar, "kA", 0.98, 0.01, 0.7, 1.2, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-4, 6.0e-3, 8.5e-3, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "kC", 1.0, 0.01, 0.9, 1.30, ierrflag); iPar++;
    junoNLMinuit->mnparm(iPar, "errGamma", 0, 0.01, 0, 0, ierrflag); iPar++;
    
    //junoNLMinuit->FixParameter(0);
    //junoNLMinuit->FixParameter(1);
    //junoNLMinuit->FixParameter(2);
    //junoNLMinuit->FixParameter(3);

    // Minimization strategy
    junoNLMinuit->SetErrorDef(1);
    arglist[0]=2;
    junoNLMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 5000; //maxCalls
    arglist[1] = 0.01; // tolerance
    junoNLMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    junoNLMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

	for(int i=0; i<m_nParameter; i++)
	{
	    junoNLMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
		//cout<<"curvalue: "<<curvalue<<"	curerror: "<<curerror<<endl;
	}

    m_DoFit = true;

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;
    delete junoNLMinuit;
    return min;
}


void junoNLChiFunction::Plot()
{
    if(!m_DoFit) {
        cout << " >>> Fitting has not been finished !<<< " << endl; return;
    }
    
    if(m_gridSearch) {
        int nPoints = 3;
        double par[nPoints];
        par[0] = final_kA; par[1] = final_kB; par[2] = final_kC;
        SetParameters(par);
    } else {
        int nPoints = m_nParameter;
        double par[nPoints];
        for(int iPoint=0; iPoint<nPoints; iPoint++) {
            par[iPoint] = m_bestFit[iPoint];
        }
        SetParameters(par);
    }

    if(junoParameters::fitGammaSources) {
        gammaNLExperiment::Plot();
    }
    if(junoParameters::fitB12) {
        m_b12Data->Plot();
    }
    if(junoParameters::fitC11) {
        m_c11Data->Plot(); // 
    }

    if(junoParameters::fitC10) {
        m_c10Data->Plot(); // 
    }

}


bool junoNLChiFunction::GridSearch()
{
    double m_like = 1000;
    // begin with some initial values
    double tmp_init_kA = 0.9;
    double tmp_init_kB = 6.0e-3;
    double tmp_init_kC = 1.;

    int tag_kA = 4;
    int tag_kB = 4;
    int tag_kC = 4;

    // parameters limits:
    double kAHigh = 1.2;
    double kALow  = 0.8;
    double kBHigh = 8.5e-3;
    double kBLow  = 5.5e-3;
    double kCHigh = 1.3;
    double kCLow  = 0.8; 

    double delta = 1e10;

    for(int init_val = 0;init_val < 1;init_val++){
        double step_kA = 0.01;
        double step_kB = 1e-4;
        double step_kC = 0.01;
        for(int iteration = 0; iteration < 100; iteration++){
            for(int bin_kA = -1; bin_kA<2; bin_kA++){
                for(int bin_kB=1; bin_kB<2; bin_kB++) {
                    for(int bin_kC=1; bin_kC<2; bin_kC++) {
                        double tmp_kA = tmp_init_kA + bin_kA*step_kA;
                        double tmp_kB = tmp_init_kB + bin_kB*step_kB;
                        double tmp_kC = tmp_init_kC + bin_kC*step_kC;
                        if(tmp_kA<=kAHigh and tmp_kA>=kALow and tmp_kB<=kBHigh and tmp_kB >= kBLow and tmp_kC<=kCHigh and tmp_kC>=kCLow) {
                            double d;
                            double par[3] = {tmp_kA, tmp_kB, tmp_kC};
                            SetParameters(par);
                            d = GetChi2();

                            if(d<delta){
                                tag_kA = bin_kA;
                                tag_kB = bin_kB;
                                tag_kC = bin_kC;
                                delta = d;
                            }
                        } else continue;
                    }
                }
            }
            tmp_init_kA = tmp_init_kA + static_cast<double>(tag_kA*step_kA);
            tmp_init_kB = tmp_init_kB + static_cast<double>(tag_kB*step_kB);
            tmp_init_kC = tmp_init_kC + static_cast<double>(tag_kC*step_kC);
            if(tag_kA==0 and tag_kB==0 and tag_kC==0) {
                step_kA /= 2.;
                step_kB /= 2.;
                step_kC /= 2.;
            }
            tag_kA = tag_kB = tag_kC = 0;
            cout << tmp_init_kA << " " << tmp_init_kB << " " << tmp_init_kC << " " << delta << endl;
            cout << "iteration: " << iteration << " " << delta << endl; 
        }
    }
    if (delta<m_like) {
        final_kA = tmp_init_kA;
        final_kB = tmp_init_kB;
        final_kC = tmp_init_kC;
    }
    m_DoFit = true;
    return true;
}


bool junoNLChiFunction::ScanContour()
{
    ofstream outfile;
    outfile.open("./chi2Files/Gamma.txt");

    double m_chi2Scan = 0;
    double m_kAmin    = final_kA; 
    double m_kBmin    = final_kB; 
    double m_kCmin    = final_kC; 

    int step_num = 20;
    double kAstep = 2*m_kALimit/step_num;
    double kBstep = 2*m_kBLimit/step_num;
    double kCstep = 2*m_kCLimit/step_num;

    for(int iStep = 0; iStep<step_num; iStep++) {
        m_kBmin = final_kB - m_kBLimit+kBstep*iStep; 
        for(int jStep=0; jStep<step_num; jStep++) {
            m_kCmin = final_kC - m_kCLimit+kCstep*jStep;
            double par[3] = {m_kAmin, m_kBmin, m_kCmin};
            SetParameters(par);
            m_chi2Scan = GetChi2();
            outfile << m_kBmin << " " << m_kCmin << " " << m_chi2Scan << endl;
        }
    }
    outfile.close();

}
