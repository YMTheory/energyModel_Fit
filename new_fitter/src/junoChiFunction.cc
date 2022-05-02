#include "junoChiFunction.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"
#include "junoParameters.hh"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

gammaResponse* junoChiFunction::gamma_array[10];
string junoChiFunction::gamma_name[10];
junoB12_simplified* junoChiFunction::b12data;
michelElectron* junoChiFunction::michel;

double junoChiFunction::m_chi2    = 0.;
double junoChiFunction::m_chi2Min = 100000;
double junoChiFunction::m_chi2B12 = 0;
double junoChiFunction::m_chi2Gam = 0;

int junoChiFunction::m_nParameter = 3;
double junoChiFunction::m_bestFit[20] = {0.};
double junoChiFunction::m_bestFitError[20] = {0.};

double junoChiFunction::gamma_amp[10];
//double junoChiFunction::init_amp [10] = {1800, 1600, 1500, 1200, 1800, 1800, 1800, 1800, 2100, 1800};    // DYB
double junoChiFunction::init_amp [10] = {1800, 1600, 1700, 1600, 1600, 1600, 1600, 1600, 1600, 1800};    //  TAo
//double junoChiFunction::init_amp [10] = {2500, 2200, 2120, 2000, 2000, 1500, 1550, 1360, 1600, 1800}; // JUNO

bool junoChiFunction::m_DoFit = false;
int junoChiFunction::m_nData;

bool junoChiFunction::m_doGamFit = junoParameters::fitGammaSources;
bool junoChiFunction::m_doB12Fit = junoParameters::fitB12;
bool junoChiFunction::m_doMichel = junoParameters::fitMichelSpec;

double junoChiFunction::m_micNPE = 7.76588e+04; 
double junoChiFunction::m_micSPE = 8.19608e+02;
double junoChiFunction::m_micNPEerr = 5.52924e+01; 
double junoChiFunction::m_micSPEerr = 2.13277e+02;
double junoChiFunction::m_resp1 = 0;
double junoChiFunction::m_resp2 = 0;

junoChiFunction::junoChiFunction() {
    
    m_nData = 0;

    if (m_doGamFit) {
        if (junoParameters::expConfig == "JUNO") {
        // JUNO normal
        double init_amp1[10] = {2500, 2200, 2120, 2000, 2000, 1500, 1550, 1360, 1600, 1800}; // JUNO
        gamma_array[m_nData] = new gammaResponse("Cs137", 100, 600, 1000);          gamma_name[m_nData] = "Cs137";    gamma_amp[m_nData] = init_amp1[0]; m_nData++;    
        gamma_array[m_nData] = new gammaResponse("Mn54", 100, 900, 1300);           gamma_name[m_nData] = "Mn54";     gamma_amp[m_nData] = init_amp1[1]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Ge68", 100, 1100, 1500);          gamma_name[m_nData] = "Ge68";     gamma_amp[m_nData] = init_amp1[2]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("K40", 100, 1750, 2250);           gamma_name[m_nData] = "K40";      gamma_amp[m_nData] = init_amp1[3]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nH", 100, 2800, 3500);            gamma_name[m_nData] = "nH";       gamma_amp[m_nData] = init_amp1[4]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Co60", 100, 3200, 3700);          gamma_name[m_nData] = "Co60";     gamma_amp[m_nData] = init_amp1[5]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmBe", 100, 6000, 6800);          gamma_name[m_nData] = "AmBe";     gamma_amp[m_nData] = init_amp1[6]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nC12", 100, 6700, 7600);          gamma_name[m_nData] = "nC12";     gamma_amp[m_nData] = init_amp1[7]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmC", 100, 8400, 9400);           gamma_name[m_nData] = "AmC";      gamma_amp[m_nData] = init_amp1[8]; m_nData++;

        }
       
        if (junoParameters::expConfig == "DYB") {
        double init_amp1[10] = {1800, 1600, 1500, 1200, 1800, 1800, 1800, 1800, 2100, 1800};    // DYB
        gamma_array[m_nData] = new gammaResponse("Cs137", 100, 70,  170);         gamma_name[m_nData] = "Cs137";    gamma_amp[m_nData] = init_amp1[0]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Mn54",  120, 100, 220);         gamma_name[m_nData] = "Mn54";     gamma_amp[m_nData] = init_amp1[1]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Ge68",  120, 130, 250);         gamma_name[m_nData] = "Ge68";     gamma_amp[m_nData] = init_amp1[2]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("K40",   160, 210, 370);         gamma_name[m_nData] = "K40";      gamma_amp[m_nData] = init_amp1[3]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nH",    100, 360, 560);         gamma_name[m_nData] = "nH";       gamma_amp[m_nData] = init_amp1[4]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Co60",  100, 400, 600);         gamma_name[m_nData] = "Co60";     gamma_amp[m_nData] = init_amp1[5]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmBe",  100, 780, 1080);         gamma_name[m_nData] = "AmBe";     gamma_amp[m_nData] = init_amp1[6]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nC12",  117, 850, 1200);         gamma_name[m_nData] = "nC12";     gamma_amp[m_nData] = init_amp1[7]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmC", 100, 1100, 1500);         gamma_name[m_nData] = "AmC";      gamma_amp[m_nData] = init_amp1[8]; m_nData++;
        }

        if(junoParameters::expConfig == "TAO") {
        double init_amp1[10] = {1800, 1600, 1700, 1600, 1600, 1600, 1600, 1600, 1600, 1800};    //  TAo
        gamma_array[m_nData] = new gammaResponse("Cs137", 100, 2400,  3100);        gamma_name[m_nData] = "Cs137";    gamma_amp[m_nData] = init_amp1[0]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Mn54",  100, 3100, 3800);         gamma_name[m_nData] = "Mn54";     gamma_amp[m_nData] = init_amp1[1]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Ge68",  100, 3800, 4500);         gamma_name[m_nData] = "Ge68";     gamma_amp[m_nData] = init_amp1[2]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("K40",   100, 5800, 6900);         gamma_name[m_nData] = "K40";      gamma_amp[m_nData] = init_amp1[3]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nH",    100, 9100, 10600);        gamma_name[m_nData] = "nH";       gamma_amp[m_nData] = init_amp1[4]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Co60",  100, 10100, 11500);       gamma_name[m_nData] = "Co60";     gamma_amp[m_nData] = init_amp1[5]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmBe",  100, 19000, 21500);       gamma_name[m_nData] = "AmBe";     gamma_amp[m_nData] = init_amp1[6]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nC12",  100, 21000, 24000);       gamma_name[m_nData] = "nC12";     gamma_amp[m_nData] = init_amp1[7]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmC", 100, 26500, 29500);         gamma_name[m_nData] = "AmC";      gamma_amp[m_nData] = init_amp1[8]; m_nData++;
        }
    }

    if (m_doB12Fit) {
        b12data = new junoB12_simplified(100, 0, 3000);
    }

    if (m_doMichel) {
        michel = new michelElectron();
    }
}


junoChiFunction::~junoChiFunction()
{
    if (m_doGamFit) {
        for (int i=0; i<m_nData; i++) {
            delete gamma_array[i];
        }
    }

    if (m_doB12Fit) {
        delete b12data;
    }

    if  (m_doMichel) {
        delete michel;
    }
}


void junoChiFunction::LoadData()
{
    if (m_doGamFit) {
        for(int i=0; i<m_nData; i++) {
            gamma_array[i]->LoadData();
        }
    }

    if (m_doB12Fit) {
        b12data->Initialize();
    }

    if (m_doMichel) {
        michel->Initialize();
    }

}


double junoChiFunction::GetChi2Michel()
{
    double chi2 = 0;
    double cerRatio = 0.067;
    double calcSPE;
    if (junoParameters::pesigmaMode == "kSeparate")
        calcSPE = TMath::Sqrt(m_resp1 * m_micNPE + m_resp2 * m_micNPE * m_micNPE);
    if (junoParameters::pesigmaMode == "kSeparate")
        calcSPE = TMath::Sqrt(m_micNPE * (1-cerRatio) + electronResponse::fElecResol->Eval(m_micNPE*cerRatio) + electronResponse::fNtotCov->Eval(m_micNPE));
    chi2 = (m_micSPE - calcSPE) * (m_micSPE - calcSPE) / m_micSPEerr / m_micSPEerr;
    return chi2;
}


double junoChiFunction::GetChi2( double maxChi2)
{
    double chi2 = 0;

    if (junoParameters::doResFit) {
        if (m_doGamFit) {
            for (int i=0; i<m_nData; i++) {
                chi2 += gamma_array[i]->GetChi2();
            }
            cout << "gamma chi2 " << chi2 << " " ;
        }

        if (m_doB12Fit) {
            double tmp = b12data->GetChi2();
            chi2 += tmp;
            cout << "B12 chi2 " << tmp << " " ;
        }

        //chi2 += GetChi2Michel();
        if (m_doMichel) {
            double tmp = michel->GetChi2();
            chi2 += tmp;
            cout << "michel chi2 " << tmp << " " ;
        }
    }

    else {
        if (m_doGamFit) {
            for (int i=0; i<m_nData; i++) {
                chi2 += gamma_array[i]->GetChi2_onlyNonl();
            }
            cout << "gamma chi2 " << chi2 << " " ;
        }
        if (m_doB12Fit) {
            double tmp = b12data->GetChi2();
            chi2 += tmp;
            cout << "B12 chi2 " << tmp << " " ;
        }
    }


    cout << "current chi2 = " << chi2 << endl;

    return chi2;
}


void junoChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


void junoChiFunction::SetParameters(double *par)
{
    int iPar = 0;

    if (junoParameters::scintillatorParameterization == kSimulation or junoParameters::scintillatorParameterization == kIntegralCalc ) {
        //electronQuench::setkA(par[iPar]);                           iPar++;
        cout << par[iPar] << " " ;
        electronQuench::setEnergyScale(par[iPar]);                  iPar++;
        cout << par[iPar] << " " ;
        electronQuench::setBirk1(par[iPar]);                        iPar++;
    }

    if (junoParameters::cerenkovMode == "kSimulationCer") {
        cout << par[iPar] << " " ;
        electronCerenkov::setkC(par[iPar]);                         iPar++;
    }
    
    if (junoParameters::cerenkovMode == "kAnalyticalCer") {
        electronCerenkov::setA2(par[iPar]);                         iPar++;
        electronCerenkov::setA3(par[iPar]);                         iPar++;
        electronCerenkov::setA4(par[iPar]);                         iPar++;
    }

    if (junoParameters::cerenkovMode == "kAnalyticalNewCer") {
        cout << par[iPar] << " " ;
        electronCerenkov::setp0(par[iPar]);                         iPar++;
        cout << par[iPar] << " " ;
        electronCerenkov::setp1(par[iPar]);                         iPar++;
        cout << par[iPar] << " " ;
        electronCerenkov::setp2(par[iPar]);                         iPar++;
        cout << par[iPar] << " " ;
        electronCerenkov::setp3(par[iPar]);                         iPar++;
        cout << par[iPar] << " " ;
        electronCerenkov::setp4(par[iPar]);                         iPar++;
        cout << par[iPar] << " " ;
        electronCerenkov::setE0(par[iPar]);                         iPar++;
    }

    if (junoParameters::doResFit) {
        if (junoParameters::pesigmaMode == "kTotal") {
            electronResponse::setra(par[iPar]);                         iPar++;
            electronResponse::setrb(par[iPar]);                         iPar++;
            electronResponse::setrc(par[iPar]);                         iPar++;
            electronResponse::SetParameters();
        }

        if (junoParameters::pesigmaMode == "kNPE") {
            electronResponse::setma               (par[iPar]);     iPar++;
            electronResponse::setmb               (par[iPar]);     iPar++;
            electronResponse::setmc               (par[iPar]);     iPar++;
            electronResponse::SetParameters();

        }

        if (junoParameters::pesigmaMode == "kNew") {
        cout << par[iPar] << " " ;
            electronResponse::setna               (par[iPar]);     iPar++;
            cout << par[iPar] << " " ;
            electronResponse::setnb               (par[iPar]);     iPar++;
            cout << par[iPar] << " " ;
            electronResponse::setnc               (par[iPar]);     iPar++;
            cout << endl;
            //electronResponse::setna1               (par[iPar]);     iPar++;
            //electronResponse::setnc1               (par[iPar]);     iPar++;
            electronResponse::SetParameters();

        }

        if (junoParameters::pesigmaMode == "kSeparate") {
            electronResponse::setc0(par[iPar]) ;                        iPar++;
            electronResponse::setc1(par[iPar]) ;                        iPar++;
            electronResponse::setc2(par[iPar]) ;                        iPar++;
            electronResponse::setd0(par[iPar]) ;                        iPar++;
            electronResponse::setd1(par[iPar]) ;                        iPar++;
            electronResponse::setd2(par[iPar]) ;                        iPar++;
            electronResponse::SetParameters();
        }

        if (m_doGamFit) {
            for (int i=0; i<m_nData; i++) {
                gamma_array[i]->SetAmp(par[iPar]);                      iPar++;
            }
        }
    }

    else {
        junoParameters::pesigmaMode = "kNew" ;
        electronResponse::setna               (par[iPar]);     iPar++;
        electronResponse::setnb               (par[iPar]);     iPar++;
        electronResponse::setnc               (par[iPar]);     iPar++;
        //electronResponse::setna1               (par[iPar]);     iPar++;
        //electronResponse::setnc1               (par[iPar]);     iPar++;
        electronResponse::SetParameters();
    }


}


double junoChiFunction::GetChiSquare(double maxChi2)
{
    junoMinuit = new TMinuit();
    junoMinuit->SetFCN(ChisqFCN);
    junoMinuit->SetPrintLevel(1);

    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    junoMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    if (junoParameters::scintillatorParameterization == kSimulation ) {
        
        //junoMinuit->mnparm(iPar, "kA", 1.00, 0.001, 0.9, 1.1, ierrflag);          iPar++;
        //junoMinuit->mnparm(iPar, "scale", 1408, 1, 1200, 1600, ierrflag);         iPar++; // JUNO Normal
        if (junoParameters::expConfig == "JUNO"){
        junoMinuit->mnparm(iPar, "scale",1406, 1, 1350, 1450, ierrflag);         iPar++; // JUNO LY / 8
        junoMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.0e-3, 7.5e-3, ierrflag);     iPar++;
        }
        
        if (junoParameters::expConfig == "DYB"){
        junoMinuit->mnparm(iPar, "scale", 205, 1, 200, 210, ierrflag);         iPar++; // JUNO LY / 8
        junoMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.0e-3, 7.5e-3, ierrflag);     iPar++;
        }
        
        if (junoParameters::expConfig == "TAO"){
        junoMinuit->mnparm(iPar, "scale", 3656, 1, 3500, 4600, ierrflag);         iPar++; // JUNO Normal
        junoMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.5e-3, 7.5e-3, ierrflag);     iPar++;
        }
    }

    if (junoParameters::scintillatorParameterization == kIntegralCalc) {
        
        //junoMinuit->mnparm(iPar, "kA", 1.00, 0.001, 0.9, 1.1, ierrflag);          iPar++;
        junoMinuit->mnparm(iPar, "scale", 1415, 1, 1200, 1600, ierrflag);         iPar++;
        //junoMinuit->FixParameter(iPar-1);
        junoMinuit->mnparm(iPar, "kB", 5.2e-3, 1e-5, 5.1e-3, 7.5e-3, ierrflag);     iPar++;
    }

    if (junoParameters::cerenkovMode == "kSimulationCer") {
        if (junoParameters::expConfig == "DYB") {
            junoMinuit->mnparm(iPar, "kC", 0.125, 0.001, 0.0, 1.0, ierrflag);            iPar++;
        }
        if (junoParameters::expConfig == "JUNO"){
            junoMinuit->mnparm(iPar, "kC", 1.00, 0.001, 0.5, 1.5, ierrflag);            iPar++;
        }
        if (junoParameters::expConfig == "TAO") {
            junoMinuit->mnparm(iPar, "kC", 2.70, 0.001, 2.0, 4.5, ierrflag);            iPar++;
        }
    }

    if (junoParameters::cerenkovMode == "kAnalyticalCer") {
        junoMinuit->mnparm(iPar, "A2", -7.90, 0.001, -20, 0, ierrflag);           iPar++;
        junoMinuit->mnparm(iPar, "A3", 13.84, 0.01, 0.0, 50, ierrflag);           iPar++;
        junoMinuit->mnparm(iPar, "A4", 0.0364, 0.0001, 0.0, 0.1, ierrflag);       iPar++;
    }

    if (junoParameters::cerenkovMode == "kAnalyticalNewCer") {
        //junoMinuit->mnparm(iPar, "p0", -0.41,   0.001,  -10., 50,  ierrflag);       iPar++;
        //junoMinuit->mnparm(iPar, "p1",  0.42,   0.001,  -10., 50,  ierrflag);       iPar++;
        //junoMinuit->mnparm(iPar, "p2", -0.016,  0.0001, -10., 50,  ierrflag);       iPar++;
        //junoMinuit->mnparm(iPar, "p3", 3.13,    0.001,  -10., 50,  ierrflag);       iPar++;
        //junoMinuit->mnparm(iPar, "p4",  0.022,  0.001,  -10., 50,  ierrflag);       iPar++;
        junoMinuit->mnparm(iPar, "p0",  0.4,    0.001,  0.2, 0.6,  ierrflag);       iPar++;
        junoMinuit->FixParameter(iPar-1);
        junoMinuit->mnparm(iPar, "p1",  0.122,   0.001, 0, 0.2,  ierrflag);       iPar++;
        junoMinuit->FixParameter(iPar-1);
        junoMinuit->mnparm(iPar, "p2",  0.573,  0.0001, 0.3, 0.8,  ierrflag);       iPar++;
        junoMinuit->FixParameter(iPar-1);
        junoMinuit->mnparm(iPar, "p3",  7.0,   0.01,    5., 10, ierrflag);       iPar++;  // JUNO LY / 8
        junoMinuit->FixParameter(iPar-1);
        //junoMinuit->mnparm(iPar, "p3",  77.2,   0.001,    0., 800, ierrflag);       iPar++;  // JUNO normal
        junoMinuit->mnparm(iPar, "p4",  0.179,  0.001,  0.1, 0.3,  ierrflag);       iPar++;
        junoMinuit->FixParameter(iPar-1);
        junoMinuit->mnparm(iPar, "E0",  0.200,  0.001,  0., 0.4,  ierrflag);          iPar++;
        junoMinuit->FixParameter(iPar-1);
    
    }

    if (junoParameters::doResFit) {

        if (junoParameters::pesigmaMode == "kTotal") {
            junoMinuit->mnparm(iPar, "ra", 0, 0.01, -20, 20, ierrflag );              iPar++;
            junoMinuit->mnparm(iPar, "rb", 1315, 1, 1250, 1850, ierrflag);            iPar++;
            junoMinuit->mnparm(iPar, "rc", 160, 1, 100, 220, ierrflag);               iPar++;
        }

        if (junoParameters::pesigmaMode == "kNPE") {
            junoMinuit->mnparm(iPar, "a", 1.02, 0.001, 0, 10, ierrflag);             iPar++;
            //junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "b", 7.9e-03, 1e-5, 1e-3, 1e-2, ierrflag);             iPar++;
            //junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "c", 0, 0.01, -1, 1, ierrflag);         iPar++;

        }

        if (junoParameters::pesigmaMode == "kNew") {
            if(junoParameters::expConfig == "DYB") {
            junoMinuit->mnparm(iPar, "a", 0.987, 1e-4, 0.8, 1.2, ierrflag);             iPar++;
            junoMinuit->mnparm(iPar, "b", 0.036, 1e-5, 0.005, 0.055, ierrflag);             iPar++;
            junoMinuit->mnparm(iPar, "n", 1.62, 0.01, 1.3, 2.2, ierrflag);          iPar++;
            }
            if(junoParameters::expConfig == "JUNO") {
            junoMinuit->mnparm(iPar, "a", 0.96, 1e-4, 0.8, 1.2, ierrflag);             iPar++;
            junoMinuit->mnparm(iPar, "b", 0.056, 1e-5, 0.005, 0.095, ierrflag);             iPar++;
            junoMinuit->mnparm(iPar, "n", 1.57, 0.01, 1.3, 2.2, ierrflag);          iPar++;
            }
            if(junoParameters::expConfig == "TAO") {
            junoMinuit->mnparm(iPar, "a", 0.803, 1e-4, 0.8, 1.5, ierrflag);             iPar++;
            //junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "b", 0.055, 1e-5, 0.001, 0.055, ierrflag);             iPar++;
            //junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "n", 1.66, 0.01, 1.3, 2.2, ierrflag);          iPar++;
            //junoMinuit->FixParameter(iPar-1);
            }
        
        }

        if (junoParameters::pesigmaMode == "kSeparate") {
            junoMinuit->mnparm(iPar, "c0", 0, 1, -200, 200, ierrflag);         iPar++;
            junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "c1", 1.77, 1e-4, 0, 50, ierrflag);            iPar++;
            junoMinuit->mnparm(iPar, "c2", 5.77e-2, 1e-4, 0, 10, ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "d0", 0, 0.01, -100, 100, ierrflag);             iPar++;
            junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "d1", 0.20708522, 1e-4, 0, 10, ierrflag);     iPar++;
            //junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "d2", 0.00237574, 1e-6, 0, 10, ierrflag);  iPar++;
            //junoMinuit->FixParameter(iPar-1);
        }

        if (m_doGamFit and junoParameters::expConfig!="TAO") {
            for(int j=0; j<m_nData; j++) {
                junoMinuit->mnparm(iPar, gamma_name[j].c_str(), gamma_amp[j], 1, gamma_amp[j]-500, gamma_amp[j]+500, ierrflag); iPar++;
            }
        }
    }

    else {
            junoParameters::pesigmaMode = "kNew";
            junoMinuit->mnparm(iPar, "a", 0.988, 1e-4, 0, 10, ierrflag);             iPar++;
            junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "b", 0.036, 1e-4, 0, 10, ierrflag);             iPar++;
            junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "n", 1.625, 0.01, 1.0, 2.0, ierrflag);          iPar++;
            junoMinuit->FixParameter(iPar-1);
    
    }

    // Minimization strategy
    junoMinuit->SetErrorDef(1);
    arglist[0]=2;
    junoMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 50000; //maxCalls
    arglist[1] = 0.01; // tolerance
    junoMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    junoMinuit->fCstatu.Data();

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    junoMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    m_nParameter = junoMinuit->GetNumPars();
	for(int i=0; i<m_nParameter; i++)
	{
	    junoMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
	}

    m_DoFit = true;

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << " with nData = " << m_nData << " and nPar = " << m_nParameter << endl;
    cout << " ====================== " << endl;

    for(int i=0; i<9; i++) {
        cout << m_bestFit[i]  << " " ;
    }
    cout << endl;
    for(int i=0; i<9; i++) {
        cout << m_bestFitError[i]  << " " ;
    }
    cout << endl;
    
    delete junoMinuit;
    return min;
}


void junoChiFunction::Plot()
{
    double par[3];
    par[0] = m_bestFit[0];
    par[1] = m_bestFit[1];
    par[2] = m_bestFit[2];

    //electronResponse::SaveElecNonl(par);

    if (m_doGamFit) {
        for (int i=0; i<m_nData; i++) {
            if(junoParameters::doResFit)
                gamma_array[i]->SaveHist();
            else
                PlotGamNonl();
        }
    }

    if (m_doB12Fit) {
        b12data->Plot();
    }


    if (m_doMichel) {
        michel->Plot();
    }
}



void junoChiFunction::PlotGamNonl()
{
    TGraphErrors* gData = new TGraphErrors();
    gData->SetName("data");
    TGraphErrors* gCalc = new TGraphErrors();
    gCalc->SetName("calc");
    
    for (int i=0; i<m_nData; i++) {
        double Etrue = gamma_array[i]->GetEtrue();
        double nonlData = gamma_array[i]->GetNonlData();
        double nonlErrData = gamma_array[i]->GetNonlErr();
        double nonlCalc = gamma_array[i]->GetNonlCalc();

        gData->SetPoint(i, Etrue, nonlData);
        gData->SetPointError(i, 0, nonlErrData);

        gCalc->SetPoint(i, Etrue, nonlCalc);
    }

    TFile* ff = new TFile("GammaNonlinearity.root", "recreate");
    gData->Write();
    gCalc->Write();
    ff->Close();

}



