#include "junoChiFunction.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"
#include "junoParameters.hh"
#include "junoSpectrum.hh"

#include <TFile.h>
#include <TGraphErrors.h>

using namespace std;

gammaResponse* junoChiFunction::Cs137;
gammaResponse* junoChiFunction::Mn54;
gammaResponse* junoChiFunction::K40;
gammaResponse* junoChiFunction::Ge68;
gammaResponse* junoChiFunction::Co60;
gammaResponse* junoChiFunction::nH;
gammaResponse* junoChiFunction::AmBe;
gammaResponse* junoChiFunction::nC12;
gammaResponse* junoChiFunction::AmC;

std::vector<gammaResponse*> junoChiFunction::gammaResponse_vector;

junoB12_simplified* junoChiFunction::b12data;

double junoChiFunction::m_chi2    = 0.;
double junoChiFunction::m_chi2Min = 100000;
double junoChiFunction::m_chi2B12 = 0;
double junoChiFunction::m_chi2Gam = 0;

int junoChiFunction::m_nParameter = 3;
double junoChiFunction::m_bestFit[20] = {0.};
double junoChiFunction::m_bestFitError[20] = {0.};

bool junoChiFunction::m_DoFit = false;

int junoChiFunction::m_nGam;
int junoChiFunction::m_nData;
std::string junoChiFunction::source_name[20];


bool junoChiFunction::m_doGamFit = junoParameters::fitGammaSources;
bool junoChiFunction::m_doB12Fit = junoParameters::fitB12;


    junoChiFunction::junoChiFunction() {

        cout << " >>>>>>>>>>>> Fitting Mode ==> " << junoParameters::m_calcOption << endl;

        m_nData = 0;
        m_nGam  = 0;

        if (m_doGamFit) {

            Cs137 = new gammaResponse("Cs137", 100, 600, 1000);
            source_name[m_nData] = "Cs137";
            gammaResponse_vector.push_back(Cs137);
            m_nData++;
            m_nGam++;

            Mn54 = new gammaResponse("Mn54", 100, 900, 1300);
            source_name[m_nData] = "Mn54";
            gammaResponse_vector.push_back(Mn54);
            m_nData++;
            m_nGam++;

            Ge68 = new gammaResponse("Ge68", 100, 1100, 1500);
            source_name[m_nData] = "Ge68";
            gammaResponse_vector.push_back(Ge68);
            m_nData++;
            m_nGam++;

            K40 = new gammaResponse("K40", 100, 1750, 2250);
            source_name[m_nData] = "K40";
            gammaResponse_vector.push_back(K40);
            m_nData++;
            m_nGam++;

            nH = new gammaResponse("nH", 100, 2800, 3500);
            source_name[m_nData] = "nH";
            gammaResponse_vector.push_back(nH);
            m_nData++;
            m_nGam++;

            Co60 = new gammaResponse("Co60", 100, 3200, 3700);
            source_name[m_nData] = "Co60";
            gammaResponse_vector.push_back(Co60);
            m_nData++;
            m_nGam++;

            AmBe = new gammaResponse("AmBe", 100, 6000, 6800);
            source_name[m_nData] = "AmBe";
            gammaResponse_vector.push_back(AmBe);
            m_nData++;
            m_nGam++;

            nC12 = new gammaResponse("nC12", 100, 6700, 7600);
            source_name[m_nData] = "nC12";
            gammaResponse_vector.push_back(nC12);
            m_nData++;
            m_nGam++;

            AmC = new gammaResponse("AmC", 100, 8400, 9400);
            source_name[m_nData] = "AmC";
            gammaResponse_vector.push_back(AmC);
            m_nData++;
            m_nGam++;


        }


        if (m_doB12Fit) {
            b12data = new junoB12_simplified(100, 4500, 17500);
        }

        electronResponse::FuncConstruct();
        electronResponse::loadElecResol();
    }

    junoChiFunction::~junoChiFunction() {
        if (m_doGamFit) {
            delete Cs137;
            delete Mn54;
            delete Ge68;
            delete K40;
            delete Co60;
            delete nH;
            delete AmBe;
            delete nC12;
            delete AmC;
        }
        if(m_doB12Fit)
            delete b12data;
    }

    void junoChiFunction::LoadData()
    {
        if (m_doGamFit) {
            for (int i=0; i<m_nGam; i++) {
                gammaResponse_vector[i]->LoadData();
            }
        }
        if (m_doB12Fit) {
            b12data->Initialize();
        }
    }




    double junoChiFunction::GetChi2( double maxChi2 )
    {
        double chi2 = 0;

        if (m_doGamFit ) {
            for(int iSource=0; iSource<m_nData; iSource++) {
                chi2 += gammaResponse_vector[iSource]->GetChi2();
            }
        }
        if(m_doB12Fit)
            chi2 += b12data->GetChi2();    

        cout << "current total chi2 = " << chi2 << endl;
        return chi2;
    }

    void junoChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
    {
        SetParameters(par);
        fval = GetChi2();
    }

    void junoChiFunction::SetParameters(double *par)
    {

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

        double min, edm, errdef;
        delete junoMinuit;
        return min;
    }



    void junoChiFunction::Plot()
    {
        if (m_doGamFit) {
        Cs137->SaveHist();
        Mn54->SaveHist();
        Ge68->SaveHist();
        K40->SaveHist();
        nH->SaveHist();
        Co60->SaveHist();
        AmBe->SaveHist();
        nC12->SaveHist();
        AmC->SaveHist();
    }
    if(m_doB12Fit)
        b12data->Plot();
}



