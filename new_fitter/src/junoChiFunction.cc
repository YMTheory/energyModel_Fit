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
        if (junoParameters::expConfig == "JUNO" ) {
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
       
        if (junoParameters::expConfig == "DYB" ) {
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

        if (junoParameters::expConfig == "Det1" ) {
        // JUNO normal
        double init_amp1[10] = {2000, 2200, 2000, 2000, 2000, 1800, 1800, 1600, 1800, 1800}; // Det1
        gamma_array[m_nData] = new gammaResponse("Cs137", 100, 700, 1000);          gamma_name[m_nData] = "Cs137";    gamma_amp[m_nData] = init_amp1[0]; m_nData++;    
        gamma_array[m_nData] = new gammaResponse("Mn54", 100, 900, 1300);           gamma_name[m_nData] = "Mn54";     gamma_amp[m_nData] = init_amp1[1]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Ge68", 100, 1100, 1500);          gamma_name[m_nData] = "Ge68";     gamma_amp[m_nData] = init_amp1[2]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("K40", 100, 1750, 2250);           gamma_name[m_nData] = "K40";      gamma_amp[m_nData] = init_amp1[3]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nH", 100, 2800, 3500);            gamma_name[m_nData] = "nH";       gamma_amp[m_nData] = init_amp1[4]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Co60", 100, 3100, 3700);          gamma_name[m_nData] = "Co60";     gamma_amp[m_nData] = init_amp1[5]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmBe", 100, 5900, 6800);          gamma_name[m_nData] = "AmBe";     gamma_amp[m_nData] = init_amp1[6]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nC12", 100, 6500, 7600);          gamma_name[m_nData] = "nC12";     gamma_amp[m_nData] = init_amp1[7]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmC", 100, 8300, 9500);           gamma_name[m_nData] = "AmC";      gamma_amp[m_nData] = init_amp1[8]; m_nData++;

        }

        if (junoParameters::expConfig == "Det3" ) {
        // JUNO normal
        double init_amp1[10] = {1800, 2000, 2000, 1700, 1700, 1600, 1700, 180, 180, 180}; // Det3
        gamma_array[m_nData] = new gammaResponse("Cs137", 100, 1300, 1800);          gamma_name[m_nData] = "Cs137";    gamma_amp[m_nData] = init_amp1[0]; m_nData++;    
        gamma_array[m_nData] = new gammaResponse("Mn54", 100, 1700, 2400);           gamma_name[m_nData] = "Mn54";     gamma_amp[m_nData] = init_amp1[1]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Ge68", 100, 2000, 2700);          gamma_name[m_nData] = "Ge68";     gamma_amp[m_nData] = init_amp1[2]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("K40", 100, 3300, 4200);           gamma_name[m_nData] = "K40";      gamma_amp[m_nData] = init_amp1[3]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nH", 100, 5350, 6450);            gamma_name[m_nData] = "nH";       gamma_amp[m_nData] = init_amp1[4]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Co60", 100, 5800, 6900);          gamma_name[m_nData] = "Co60";     gamma_amp[m_nData] = init_amp1[5]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmBe", 100, 11200, 13200);          gamma_name[m_nData] = "AmBe";     gamma_amp[m_nData] = init_amp1[6]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nC12", 100, 13000, 15000);          gamma_name[m_nData] = "nC12";     gamma_amp[m_nData] = init_amp1[7]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmC", 100, 16800, 18800);           gamma_name[m_nData] = "AmC";      gamma_amp[m_nData] = init_amp1[8]; m_nData++;
        }

        if (junoParameters::expConfig == "Det5" ) {
        // JUNO normal
        double init_amp1[10] = {2100, 1800, 1600, 1300, 1800, 1900, 1600, 1400, 1800, 1800}; // Det3
        gamma_array[m_nData] = new gammaResponse("Cs137", 80, 70, 150);             gamma_name[m_nData] = "Cs137";    gamma_amp[m_nData] = init_amp1[0]; m_nData++;    
        gamma_array[m_nData] = new gammaResponse("Mn54", 100, 90, 190);          gamma_name[m_nData] = "Mn54";     gamma_amp[m_nData] = init_amp1[1]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Ge68", 100, 110, 210);          gamma_name[m_nData] = "Ge68";     gamma_amp[m_nData] = init_amp1[2]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("K40", 130, 180, 310);           gamma_name[m_nData] = "K40";      gamma_amp[m_nData] = init_amp1[3]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nH", 100, 300, 500);            gamma_name[m_nData] = "nH";       gamma_amp[m_nData] = init_amp1[4]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Co60", 90, 330, 510);          gamma_name[m_nData] = "Co60";     gamma_amp[m_nData] = init_amp1[5]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmBe", 120, 660, 900);          gamma_name[m_nData] = "AmBe";     gamma_amp[m_nData] = init_amp1[6]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nC12", 115, 760, 990);          gamma_name[m_nData] = "nC12";     gamma_amp[m_nData] = init_amp1[7]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmC", 100, 940, 1240);           gamma_name[m_nData] = "AmC";      gamma_amp[m_nData] = init_amp1[8]; m_nData++;

        }

        if (junoParameters::expConfig == "Det6" ) {
        // JUNO normal
        double init_amp1[10] = {1800, 1700, 1500, 1700, 1800, 1800, 1800, 1800, 2000, 2000}; // Det3
        gamma_array[m_nData] = new gammaResponse("Cs137", 100, 330, 530);             gamma_name[m_nData] = "Cs137";    gamma_amp[m_nData] = init_amp1[0]; m_nData++;    
        gamma_array[m_nData] = new gammaResponse("Mn54", 100, 450, 650);          gamma_name[m_nData] = "Mn54";     gamma_amp[m_nData] = init_amp1[1]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Ge68", 100, 550, 750);          gamma_name[m_nData] = "Ge68";     gamma_amp[m_nData] = init_amp1[2]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("K40", 100, 850, 1150);           gamma_name[m_nData] = "K40";      gamma_amp[m_nData] = init_amp1[3]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nH", 100, 1350, 1750);            gamma_name[m_nData] = "nH";       gamma_amp[m_nData] = init_amp1[4]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Co60", 100, 1500, 1900);          gamma_name[m_nData] = "Co60";     gamma_amp[m_nData] = init_amp1[5]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmBe", 100, 2900, 3500);          gamma_name[m_nData] = "AmBe";     gamma_amp[m_nData] = init_amp1[6]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nC12", 100, 3200, 3900);          gamma_name[m_nData] = "nC12";     gamma_amp[m_nData] = init_amp1[7]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmC", 100, 4000, 4800);           gamma_name[m_nData] = "AmC";      gamma_amp[m_nData] = init_amp1[8]; m_nData++;

        }

        if (junoParameters::expConfig == "Det7" ) {
        // JUNO normal
        double init_amp1[10] = {1328, 2347, 2156, 1712, 1994, 1983, 1807, 1624, 1925, 1800}; // Det3
        gamma_array[m_nData] = new gammaResponse("Cs137", 150, 150, 300);             gamma_name[m_nData] = "Cs137";    gamma_amp[m_nData] = init_amp1[0]; m_nData++;    
        gamma_array[m_nData] = new gammaResponse("Mn54", 100, 180, 380);          gamma_name[m_nData] = "Mn54";     gamma_amp[m_nData] = init_amp1[1]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Ge68", 100, 230, 430);          gamma_name[m_nData] = "Ge68";     gamma_amp[m_nData] = init_amp1[2]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("K40", 100, 400, 600);           gamma_name[m_nData] = "K40";      gamma_amp[m_nData] = init_amp1[3]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nH", 100, 650, 950);            gamma_name[m_nData] = "nH";       gamma_amp[m_nData] = init_amp1[4]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("Co60", 100, 650, 950);          gamma_name[m_nData] = "Co60";     gamma_amp[m_nData] = init_amp1[5]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmBe", 100, 1400, 1800);          gamma_name[m_nData] = "AmBe";     gamma_amp[m_nData] = init_amp1[6]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("nC12", 100, 1550, 1950);          gamma_name[m_nData] = "nC12";     gamma_amp[m_nData] = init_amp1[7]; m_nData++;
        gamma_array[m_nData] = new gammaResponse("AmC", 100, 2000, 2500);           gamma_name[m_nData] = "AmC";      gamma_amp[m_nData] = init_amp1[8]; m_nData++;

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
        electronQuench::setEnergyScale(par[iPar]);                  iPar++;
        cout << par[iPar-1] << " " ;
        electronQuench::setBirk1(par[iPar]);                        iPar++;
        cout << par[iPar-1] << " " ;
    }

    if (junoParameters::cerenkovMode == "kSimulationCer" or junoParameters::cerenkovMode=="kSimulationCerLocal") {
        electronCerenkov::setkC(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
    }
    
    if (junoParameters::cerenkovMode == "kAnalyticalCer") {
        electronCerenkov::setA2(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setA3(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setA4(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setE0(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
    }

    if (junoParameters::cerenkovMode == "kAnalyticalNewCer") {
        electronCerenkov::setp0(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setp1(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setp2(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setp3(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setp4(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setE0(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
    }

    if (junoParameters::cerenkovMode == "kAnalyticalNewCer1") {
        electronCerenkov::setp0(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setp1(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setp2(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
        electronCerenkov::setE0(par[iPar]);                         iPar++;
        cout << par[iPar-1] << " " ;
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
    int nPhysPar = 0;
    junoMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    if (junoParameters::scintillatorParameterization == kSimulation ) {
        
        //junoMinuit->mnparm(iPar, "kA", 1.00, 0.001, 0.9, 1.1, ierrflag);          iPar++;
        //junoMinuit->mnparm(iPar, "scale", 1408, 1, 1200, 1600, ierrflag);         iPar++; // JUNO Normal
        if (junoParameters::expConfig == "JUNO" or junoParameters::expConfig=="Det1"){
        junoMinuit->mnparm(iPar, "scale",1400, 1, 1350, 1450, ierrflag);         iPar++; // JUNO LY / 8
        junoMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.5e-3, 7.5e-3, ierrflag);     iPar++;
        }

        if (junoParameters::expConfig=="Det3"){
        junoMinuit->mnparm(iPar, "scale",2800, 1, 2550, 3150, ierrflag);            iPar++; // JUNO LY / 8
        junoMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.0e-3, 7.5e-3, ierrflag);     iPar++;
        }
        
        if (junoParameters::expConfig=="Det6"){
        junoMinuit->mnparm(iPar, "scale",700, 1, 500, 800, ierrflag);            iPar++; // JUNO LY / 8
        junoMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.0e-3, 7.5e-3, ierrflag);     iPar++;
        }
        if (junoParameters::expConfig=="Det7"){
        junoMinuit->mnparm(iPar, "scale", 347.4, 1, 250, 450, ierrflag);            iPar++; // JUNO LY / 8
        junoMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-5, 5.0e-3, 7.5e-3, ierrflag);     iPar++;
        }
        
        if (junoParameters::expConfig == "DYB" or junoParameters::expConfig == "Det5"){
        junoMinuit->mnparm(iPar, "scale", 174.3, 1, 150, 200, ierrflag);         iPar++; // JUNO LY / 8
        junoMinuit->mnparm(iPar, "kB", 6.0e-3, 1e-5, 5.0e-3, 7.5e-3, ierrflag);     iPar++;
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

    if (junoParameters::cerenkovMode == "kSimulationCer" or junoParameters::cerenkovMode=="kSimulationCerLocal") {
        if (junoParameters::expConfig == "DYB" or junoParameters::expConfig == "Det5") {
            junoMinuit->mnparm(iPar, "kC", 0.124, 0.001, 0.0, 1.0, ierrflag);            iPar++;
        }
        if (junoParameters::expConfig == "JUNO" or junoParameters::expConfig=="Det1"){
            junoMinuit->mnparm(iPar, "kC", 1.00, 0.001, 0.5, 1.5, ierrflag);            iPar++;
        }
        if ( junoParameters::expConfig=="Det6"){
            junoMinuit->mnparm(iPar, "kC", 0.50, 0.001, 0.0, 1.5, ierrflag);            iPar++;
        }
        if ( junoParameters::expConfig=="Det7"){
            junoMinuit->mnparm(iPar, "kC", 0.265, 0.001, 0.0, 1.5, ierrflag);            iPar++;
        }
        if (junoParameters::expConfig == "TAO") {
            junoMinuit->mnparm(iPar, "kC", 2.70, 0.001, 2.0, 4.5, ierrflag);            iPar++;
        }
        if (junoParameters::expConfig == "Det3") {
            junoMinuit->mnparm(iPar, "kC", 2.00, 0.001, 2.0, 4.5, ierrflag);            iPar++;
        }
    }

    if (junoParameters::cerenkovMode == "kAnalyticalCer") {
        if (junoParameters::expConfig == "Det1") {
            junoMinuit->mnparm(iPar, "A2", -10.20, 0.001, -20, 0, ierrflag);          iPar++;
            junoMinuit->mnparm(iPar, "A3", 14.33, 0.01, 0.0, 30, ierrflag);           iPar++;
            junoMinuit->mnparm(iPar, "A4", 0.0364, 0.0001, 0.0, 0.1, ierrflag);       iPar++;
            junoMinuit->mnparm(iPar, "E0",  0.183,  0.001,  0., 0.3,  ierrflag);      iPar++;
        }
        if (junoParameters::expConfig == "Det5") {
            junoMinuit->mnparm(iPar, "A2", -1.275, 0.001, -20, 0, ierrflag);          iPar++;
            junoMinuit->mnparm(iPar, "A3", 1.79, 0.01, 0.0, 10, ierrflag);           iPar++;
            junoMinuit->mnparm(iPar, "A4", 0.0364, 0.0001, 0.0, 0.1, ierrflag);       iPar++;
            junoMinuit->mnparm(iPar, "E0",  0.183,  0.001,  0., 0.3,  ierrflag);      iPar++;
        }
    }

    if (junoParameters::cerenkovMode == "kAnalyticalNewCer") {

        // JUNO Sim
        if(junoParameters::expConfig == "JUNO") {
            junoMinuit->mnparm(iPar, "p0",  0.9,    0.001,  0.6, 1.2,  ierrflag);       iPar++;
            junoMinuit->mnparm(iPar, "p1",  10.21,   0.001, 5, 15,  ierrflag);       iPar++;   // Det1
            junoMinuit->mnparm(iPar, "p2",  0.010,  0.001, 0., 0.02,  ierrflag);       iPar++;   // Det1
            junoMinuit->mnparm(iPar, "p3",  71.6,   0.001,    50., 100, ierrflag);       iPar++;  // JUNO normal
            junoMinuit->mnparm(iPar, "p4",  -9.9,  0.001,  -15, -5,  ierrflag);       iPar++;
            junoMinuit->mnparm(iPar, "E0",  0.134,  0.001,  0., 0.3,  ierrflag);          iPar++;
        }
        
        // Det1 Sim
        if (junoParameters::expConfig == "Det1") {
            junoMinuit->mnparm(iPar, "p0", 121 , 0.1, 100 ,  140 , ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "p1", 0.68, 0.001, 0.4   ,0.9  , ierrflag);      iPar++;
            //junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "p2", -0.083 , 0.001, -0.1 , -0.04, ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "p3", 73.410, 0.001, 60  , 90  , ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "p4", -9.82 , 0.001, -15 , -5  , ierrflag);      iPar++;
            //junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "E0", 0.267 , 0.001, 0.1 , 0.4 , ierrflag);      iPar++;
            //junoMinuit->FixParameter(iPar-1);
        }

        //junoMinuit->FixParameter(iPar-1);
        //junoMinuit->mnparm(iPar, "p1",  0.122,   0.001, 0, 0.2,  ierrflag);       iPar++; 
        //junoMinuit->FixParameter(iPar-1);
        //junoMinuit->mnparm(iPar, "p2",  0.573,  0.0001, 0.3, 0.8,  ierrflag);       iPar++;
        //junoMinuit->FixParameter(iPar-1);
        //junoMinuit->mnparm(iPar, "p3",  7.0,   0.01,    5., 10, ierrflag);       iPar++;  // JUNO LY / 8
        //junoMinuit->mnparm(iPar, "p3",  49,   0.001,    0., 800, ierrflag);       iPar++;  // Det1
        //junoMinuit->FixParameter(iPar-1);
        //junoMinuit->FixParameter(iPar-1);
    
    }

    if (junoParameters::cerenkovMode == "kAnalyticalNewCer1" ) {
        if (junoParameters::expConfig == "JUNO") {
            junoMinuit->mnparm(iPar, "p0", 121 , 0.1, 100 ,  140 , ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "p1", 0.54, 0.001, 0.4   ,0.9  , ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "p2", -0.130 , 0.001, -0.4 , -0.04, ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "E0", 0.154 , 0.001, 0.1 , 0.4 , ierrflag);      iPar++;
        }
        if (junoParameters::expConfig == "Det1") {
            junoMinuit->mnparm(iPar, "p0", 121 , 0.1, 100 ,  140 , ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "p1", 0.68, 0.001, 0.4   ,0.9  , ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "p2", 0.083 , 0.001, -0.4 , 0.4, ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "E0", 0.267 , 0.001, 0.1 , 0.4 , ierrflag);      iPar++;
        }
        if (junoParameters::expConfig == "Det5") {
            junoMinuit->mnparm(iPar, "p0", 17.2 , 0.1, 0 ,  40 , ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "p1", 1.50, 0.001, 0.4   ,5.5  , ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "p2", -0.162 , 0.001, -0.4 , -0.04, ierrflag);      iPar++;
            junoMinuit->mnparm(iPar, "E0", 0.137 , 0.001, 0.1 , 0.4 , ierrflag);      iPar++;
        }
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
            if( junoParameters::expConfig == "Det6") {
            junoMinuit->mnparm(iPar, "a", 0.990, 1e-4, 0.8, 1.2, ierrflag);             iPar++;
            junoMinuit->mnparm(iPar, "b", 0.088, 1e-5, 0.005, 0.155, ierrflag);             iPar++;
            junoMinuit->mnparm(iPar, "n", 1.50, 0.01, 1.0, 2.2, ierrflag);          iPar++;
            }
            if( junoParameters::expConfig == "Det7") {
            junoMinuit->mnparm(iPar, "a", 0.989, 1e-4, 0.5, 1.5, ierrflag);             iPar++;
            junoMinuit->mnparm(iPar, "b", 0.045, 1e-5, 0.005, 0.205, ierrflag);             iPar++;
            junoMinuit->mnparm(iPar, "n", 1.500, 0.01, 1.0, 2.2, ierrflag);          iPar++;
            }
            if(junoParameters::expConfig == "DYB" or junoParameters::expConfig == "Det5") {
            junoMinuit->mnparm(iPar, "a", 0.987, 1e-4, 0.8, 1.2, ierrflag);             iPar++;
            junoMinuit->mnparm(iPar, "b", 0.092, 1e-5, 0.005, 1.0, ierrflag);             iPar++;
            junoMinuit->mnparm(iPar, "n", 1.30, 0.01, 1.0, 2.2, ierrflag);          iPar++;
            }
            if(junoParameters::expConfig == "JUNO" or junoParameters::expConfig == "Det1") {
            junoMinuit->mnparm(iPar, "a", 0.991, 1e-4, 0.8, 1.2, ierrflag);             iPar++;
            //junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "b", 0.034, 1e-5, 0.005, 0.095, ierrflag);             iPar++;
            //junoMinuit->FixParameter(iPar-1);
            junoMinuit->mnparm(iPar, "n", 1.59, 0.01, 1.2, 2.0, ierrflag);          iPar++;
            //junoMinuit->FixParameter(iPar-1);
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
        nPhysPar = iPar;

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

    for(int i=0; i<nPhysPar; i++) {
        cout << m_bestFit[i]  << " " ;
    }
    cout << endl;
    for(int i=0; i<nPhysPar; i++) {
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



