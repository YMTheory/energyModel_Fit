#include "gammaResChiFunction.hh"
#include "electronResol.hh"
#include "gammaResol.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "junoParameters.hh"

#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include <TMinuit.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH2.h>

using namespace std;

gammaResol* gammaResChiFunction::Cs137Data = 0;
gammaResol* gammaResChiFunction::Mn54Data = 0;
gammaResol* gammaResChiFunction::K40Data = 0;
gammaResol* gammaResChiFunction::nHData = 0;
gammaResol* gammaResChiFunction::Co60Data = 0;
gammaResol* gammaResChiFunction::Tl208Data = 0;
gammaResol* gammaResChiFunction::nC12Data = 0;
gammaResol* gammaResChiFunction::O16Data = 0;
gammaResol* gammaResChiFunction::nFe56Data = 0;

bool gammaResChiFunction::m_DoFit = false;
bool gammaResChiFunction::m_gridSearch = false;

double gammaResChiFunction::m_chi2    = 0.;
int gammaResChiFunction::m_nData   = 0;
double gammaResChiFunction::final_pA = 0;
double gammaResChiFunction::final_pB = 0;
double gammaResChiFunction::final_pC = 0;
double gammaResChiFunction::final_kA = 0;
double gammaResChiFunction::final_kB = 0;
double gammaResChiFunction::final_kC = 0;

int gammaResChiFunction::m_nParameter = 3;
double gammaResChiFunction::m_bestFit[20] = {0.};
double gammaResChiFunction::m_bestFitError[20] = {0.};

std::map<std::string, gammaResol*> gammaResChiFunction::mapGammaResol;
std::string gammaResChiFunction::source_name[20];


gammaResChiFunction::gammaResChiFunction() {
    Cs137Data = new gammaResol("Cs137", 700, 1100, 100);
    Cs137Data->LoadData(); 
    source_name[m_nData] = "Cs137"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Cs137", Cs137Data));

    Mn54Data  = new gammaResol("Mn54", 900, 1300, 100);
    Mn54Data->LoadData(); 
    source_name[m_nData] = "Mn54"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Mn54", Mn54Data));

    K40Data  = new gammaResol("K40", 900, 1300, 100);
    K40Data->LoadData(); 
    source_name[m_nData] = "K40"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("K40", K40Data));

    nHData  = new gammaResol("nH", 900, 1300, 100);
    nHData->LoadData(); 
    source_name[m_nData] = "nH"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("nH", nHData));

    Co60Data  = new gammaResol("Co60", 900, 1300, 100);
    Co60Data->LoadData(); 
    source_name[m_nData] = "Co60"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Co60", Co60Data));

    Tl208Data  = new gammaResol("Tl208", 900, 1300, 100);
    Tl208Data->LoadData(); 
    source_name[m_nData] = "Tl208"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Tl208", Tl208Data));

    nC12Data  = new gammaResol("nC12", 900, 1300, 100);
    nC12Data->LoadData(); 
    source_name[m_nData] = "nC12"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("nC12", nC12Data));

    O16Data  = new gammaResol("O16", 900, 1300, 100);
    O16Data->LoadData(); 
    source_name[m_nData] = "O16"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("O16", O16Data));

    nFe56Data  = new gammaResol("nFe56", 900, 1300, 100);
    nFe56Data->LoadData(); 
    source_name[m_nData] = "nFe56"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("nFe56", nFe56Data));

}

void gammaResChiFunction::Initialize()
{
    string mode;
    if   (junoParameters::gammaNLOption == 0) 
        mode = "PrimElecDist";
    else if(junoParameters::gammaNLOption == 1) 
        mode = "2-layer sampling";

    cout << "======> Nonlineairity Fitting Mode: " << mode << endl;
}



double gammaResChiFunction::GetChi2 (double maxChi2) {
    m_chi2 = 0;
    m_chi2 += Cs137Data->GetChi2();
    m_chi2 += Mn54Data->GetChi2();
    m_chi2 += K40Data->GetChi2();
    m_chi2 += nHData->GetChi2();
    m_chi2 += Co60Data->GetChi2();
    m_chi2 += Tl208Data->GetChi2();
    m_chi2 += nC12Data->GetChi2();
    m_chi2 += O16Data->GetChi2();
    m_chi2 += nFe56Data->GetChi2();

    // Add Penalty Term: 
    m_chi2 += TMath::Power(gammaResol::GetGammaScale()/junoParameters::m_gammaError, 2) ;

    //cout << "kA : " << electronQuench::getkA() << "  kB: " << electronQuench::getBirk1() << "  kC: " << electronCerenkov::getkC()  << "  chi2: " << m_chi2 << endl;
    return m_chi2;
}


void gammaResChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


void gammaResChiFunction::SetParameters(double* par)
{
    electronQuench::setkA(par[0]);
    electronQuench::setBirk1(par[1]);
    electronCerenkov::setkC(par[2]);
    gammaResol::SetGammaScale(par[3]);
    //electronResol::setpA(par[0]);
    //electronResol::setpB(par[1]);
    //electronResol::setpC(par[2]);
}

double gammaResChiFunction::GetChiSquare(double maxChi2)
{
    gammaResMinuit = new TMinuit();
    gammaResMinuit->SetFCN(ChisqFCN);
    gammaResMinuit->SetPrintLevel(1);
    
    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    gammaResMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // Configurate parameters
    gammaResMinuit->mnparm(iPar, "kA", 0.960, 0.0001, 0.95, 0.97, ierrflag); iPar++;
    gammaResMinuit->mnparm(iPar, "kB", 6.5e-3, 1e-4, 5.5e-3, 7.2e-3, ierrflag); iPar++;
    gammaResMinuit->mnparm(iPar, "kC", 1.00, 0.01, 0.98, 1.02, ierrflag); iPar++;
    gammaResMinuit->mnparm(iPar, "gammaScale", 0.003, 0.0001, 0, 0, ierrflag); iPar++;
    //gammaResMinuit->mnparm(iPar, "pA", 0.0258, 0.0001, 0.012, 0.029, ierrflag); iPar++;
    //gammaResMinuit->mnparm(iPar, "pB", 6.80e-3, 1e-5, 6.5e-3, 7.0e-3, ierrflag); iPar++;
    //gammaResMinuit->mnparm(iPar, "pC", 0.0, 0.001, 0.0, 0.00, ierrflag); iPar++;
    
    //gammaResMinuit->FixParameter(0);
    //gammaResMinuit->FixParameter(1);
    //gammaResMinuit->FixParameter(2);
    //gammaResMinuit->FixParameter(3);
    //gammaResMinuit->FixParameter(4);
    //gammaResMinuit->FixParameter(5);

    // Minimization strategy
    gammaResMinuit->SetErrorDef(1);
    arglist[0]=2;
    gammaResMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    //arglist[0] = 5000; //maxCalls
    //arglist[1] = 1.0; // tolerance
    //gammaResMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);
 
    //gammaResMinuit->mnexcm("CALL FCN", arglist, 1, ierrflag);

    arglist[0] = 5000;
    arglist[1] = 50;
    gammaResMinuit->mnexcm("MIGRAD", arglist, 2, ierrflag);



    double min, edm, errdef;
    int nvpar, nparx, icstat;
    gammaResMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

	for(int i=0; i<m_nParameter; i++)
	{
	    gammaResMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
		//cout<<"curvalue: "<<curvalue<<"	curerror: "<<curerror<<endl;
	}

    m_DoFit = true;

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;
    delete gammaResMinuit;
    return min;
}

/*bool gammaResChiFunction::GridSearch()
{
    double m_like = 1000;
    double tmp_init_pA = 0.0260;
    double tmp_init_pB = 0.0068;
    double tmp_init_pC = 0.;
    double tmp_init_kA = 0.95;
    double tmp_init_kB = 6.5e-3;
    double tmp_init_kC = 1.0;

    int tag_pA = 4;
    int tag_pB = 4;
    int tag_pC = 4;
    int tag_kA = 4;
    int tag_kB = 4;
    int tag_kC = 4;

    // parameters limits:
    double pAHigh = 0.040;
    double pALow  = 0.010;
    double pBHigh = 9.0e-3;
    double pBLow  = 5.0e-3;
    double pCHigh = 0.01;
    double pCLow  = 0.00; 
    double kAHigh = 0.90;
    double kALow  = 1.10;
    double kBHigh = 6.0e-3;
    double kBLow  = 7.0e-3;
    double kCHigh = 0.90;
    double kCLow  = 1.1; 

    double delta = 1e10;

    for(int init_val = 0;init_val < 1;init_val++){
        double step_pA = 0.001;
        double step_pB = 1e-4;
        double step_pC = 0.001;
        double step_kA = 0.01;
        double step_kB = 1e-4;
        double step_kC = 0.01;
        for(int iteration = 0; iteration < 200; iteration++){
            for(int bin_kA = -1; bin_kA<2; bin_kA++){
                for(int bin_kB=1; bin_kB<2; bin_kB++) {
                    for(int bin_kC=1; bin_kC<2; bin_kC++) {
                        double tmp_kA = tmp_init_kA + bin_kA*step_kA;
                        double tmp_kB = tmp_init_kB + bin_kB*step_kB;
                        double tmp_kC = tmp_init_kC + bin_kC*step_kC;
                        if(tmp_kA<=kAHigh and tmp_kA>=kALow and tmp_kB<=kBHigh and tmp_kB >= kBLow and tmp_kC<=kCHigh and tmp_kC>=kCLow) {
                            double d = delta;
                            double kAr[3] = {tmp_kA, tmp_kB, tmp_kC};
                            SetParameters(kAr);
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
    m_gridSearch = true;
    return true;

}
*/


void gammaResChiFunction::Plot()
{

    //if(!m_DoFit) {
    //    cout << " >>> Fitting has not been finished !<<< " << endl; return;
    //}
    //cout << " >>> Plot Nonlinearity and Resolution Fitting Results <<< " << endl;

    //if(m_gridSearch) {
    //    int nPoints = 3;
    //    double par[nPoints];
    //    par[0] = final_pA; par[1] = final_pB; par[2] = final_pC;
    //    SetParameters(par);
    //} else {
    //    int nPoints = m_nParameter;
    //    double par[nPoints];
    //    for(int iPoint=0; iPoint<nPoints; iPoint++) {
    //        par[iPoint] = m_bestFit[iPoint];
    //    }
    //    SetParameters(par);
    //}

    //double par[6] = {0.961, 6.13e-3, 1.00, 0.27, 0.0068, 0};
    //SetParameters(par);

    //cout << "Current parameters: " << electronQuench::getkA() << " " << electronQuench::getBirk1() << " " << electronCerenkov::getkC() << endl;

    TGraphErrors* gNonlData = new TGraphErrors();
    TGraph* gNonlCalc = new TGraph();
    gNonlData->SetName("gNonlData");
    gNonlCalc->SetName("gNonlCalc");

    int index = 0;
    for(int iData=0; iData<m_nData; iData++) {
        std::string source = source_name[iData];
        gammaResol* tmpGammaData = mapGammaResol.find(source)->second;
        //tmpGammaData->calcGammaNPE();
        double tmp_E       = tmpGammaData->GetEtrue();
        double tmp_pred    = tmpGammaData->GetNonlPred();
        double tmp_data    = tmpGammaData->GetNonlData();
        double tmp_dataErr = tmpGammaData->GetNonlDataErr(); 
        //cout << tmp_E << " " << tmp_data << " " << tmp_pred << endl;
        gNonlData->SetPoint(index, tmp_E, tmp_data);
        gNonlData->SetPointError(index, 0, tmp_dataErr);
        gNonlCalc->SetPoint(index, tmp_E, tmp_pred);
        index++;
    }

    gNonlData->SetMarkerStyle(20);
    gNonlData->SetMarkerColor(kBlue+1);
    gNonlData->SetLineColor(kBlue+1);
    gNonlData->SetLineWidth(3);
    gNonlData->SetMarkerSize(1.2);
    gNonlCalc->SetMarkerStyle(21);
    gNonlCalc->SetMarkerColor(kRed+1);
    gNonlCalc->SetMarkerSize(1.2);
    gNonlCalc->SetLineColor(kRed+1);
    gNonlCalc->SetLineWidth(3);

    TCanvas* c1 = new TCanvas("Nonlinearity", "Nonlinearity");
    c1->cd(); c1->SetGrid();
    gNonlData->SetTitle("Nonlinearity Fitting; Etrue/MeV; Nonlinearity");
    //gNonlData->GetYaxis()->SetRangeUser(0.01,0.045);
    gNonlData->Draw("APL");
    gNonlCalc->Draw("P SAME");
    TLegend* led = new TLegend();
    led->SetFillColor(kWhite);
    led->AddEntry(gNonlData, "data", "PL");
    led->AddEntry(gNonlCalc, "calc", "PL");
    led->Draw("SAME");


    TGraphErrors* gResData = new TGraphErrors();
    TGraph* gResCalc = new TGraph();
    gResData->SetName("gResData");
    gResCalc->SetName("gResCalc");

    index = 0;
    for(int iData=0; iData<m_nData; iData++) {
        std::string source = source_name[iData];
        gammaResol* tmpGammaData = mapGammaResol.find(source)->second;
        //tmpGammaData->calcGammaNPE();
        double tmp_E       = tmpGammaData->GetEvis();
        double tmp_pred    = tmpGammaData->GetResolPred();
        double tmp_data    = tmpGammaData->GetResolData();
        double tmp_dataErr = tmpGammaData->GetResolDataErr(); 
        //cout << tmp_E << " " << tmp_data << " " << tmp_pred << endl;
        gResData->SetPoint(index, tmp_E, tmp_data);
        gResData->SetPointError(index, 0, tmp_dataErr);
        gResCalc->SetPoint(index, tmp_E, tmp_pred);
        index++;
    }

    gResData->SetMarkerStyle(20);
    gResData->SetMarkerColor(kBlue+1);
    gResData->SetLineColor(kBlue+1);
    gResData->SetLineWidth(3);
    gResData->SetMarkerSize(1.2);
    gResCalc->SetMarkerStyle(21);
    gResCalc->SetMarkerColor(kRed+1);
    gResCalc->SetMarkerSize(1.2);
    gResCalc->SetLineColor(kRed+1);
    gResCalc->SetLineWidth(3);

    TCanvas* c2 = new TCanvas("Resolution", "Resolution");
    c2->cd(); c2->SetGrid();
    gResData->SetTitle("Resolution Fitting; Evis/MeV; Resinearity");
    //gResData->GetYaxis()->SetRangeUser(0.01,0.045);
    gResData->Draw("APL");
    gResCalc->Draw("PL SAME");
    TLegend* le = new TLegend();
    le->SetFillColor(kWhite);
    le->AddEntry(gResData, "data", "PL");
    le->AddEntry(gResCalc, "calc", "PL");
    le->Draw("SAME");



    TFile* file = new TFile("./output/gamma/gammaAllFit.root", "recreate");
    gNonlData->Write();
    gNonlCalc->Write();
    gResData->Write();
    gResCalc->Write();
    //c1->Write();
    //c2->Write();
    file->Close();
}


void gammaResChiFunction::Contour()
{
    TH2D* cont = new TH2D("cont", "", 50, 0.96, 0.97, 20, 0.9,1.10);
    for(int iA=0; iA<50; iA++) {
        double kA = 0.960+0.0002*iA;
        for(int iB=0; iB<20; iB++) {
            double kB = 6.5e-3;
            double kC = 0.9+0.01*iB;
            double par[3] = {kA, kB, kC};
            SetParameters(par);
            double chi2 = GetChi2();
            cout << kA << " " << kB << " " << kC << " " << chi2 << endl;
            cont->SetBinContent(iA, iB, chi2);
        }
    }

    TFile* out = new TFile("scan2.root", "recreate");
    cont->Write();
    out->Close();
}


