#include "gammaResChiFunction.hh"
#include "electronResol.hh"
#include "gammaResol.hh"

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

int gammaResChiFunction::m_nParameter = 3;
double gammaResChiFunction::m_bestFit[20] = {0.};
double gammaResChiFunction::m_bestFitError[20] = {0.};

std::map<std::string, gammaResol*> gammaResChiFunction::mapGammaResol;
std::string gammaResChiFunction::source_name[20];


gammaResChiFunction::gammaResChiFunction() {
    Cs137Data = new gammaResol("Cs137", 700, 1100, 100);
    Cs137Data->calcGammaNPE(); 
    source_name[m_nData] = "Cs137"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Cs137", Cs137Data));

    Mn54Data  = new gammaResol("Mn54", 900, 1300, 100);
    Mn54Data->calcGammaNPE(); 
    source_name[m_nData] = "Mn54"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Mn54", Mn54Data));

    K40Data  = new gammaResol("K40", 900, 1300, 100);
    K40Data->calcGammaNPE(); 
    source_name[m_nData] = "K40"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("K40", K40Data));

    nHData  = new gammaResol("nH", 900, 1300, 100);
    nHData->calcGammaNPE(); 
    source_name[m_nData] = "nH"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("nH", nHData));

    Co60Data  = new gammaResol("Co60", 900, 1300, 100);
    Co60Data->calcGammaNPE(); 
    source_name[m_nData] = "Co60"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("Co60", Co60Data));

    ////Tl208Data  = new gammaResol("Tl208", 900, 1300, 100);
    ////Tl208Data->calcGammaNPE(); 
    ////source_name[m_nData] = "Tl208"; m_nData++;
    ////mapGammaResol.insert(pair<std::string, gammaResol*> ("Tl208", Tl208Data));

    nC12Data  = new gammaResol("nC12", 900, 1300, 100);
    nC12Data->calcGammaNPE(); 
    source_name[m_nData] = "nC12"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("nC12", nC12Data));

    O16Data  = new gammaResol("O16", 900, 1300, 100);
    O16Data->calcGammaNPE(); 
    source_name[m_nData] = "O16"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("O16", O16Data));

    nFe56Data  = new gammaResol("nFe56", 900, 1300, 100);
    nFe56Data->calcGammaNPE(); 
    source_name[m_nData] = "nFe56"; m_nData++;
    mapGammaResol.insert(pair<std::string, gammaResol*> ("nFe56", nFe56Data));

}

double gammaResChiFunction::GetChi2 (double maxChi2) {
    m_chi2 = 0;
    m_chi2 += Cs137Data->GetChi2();
    m_chi2 += Mn54Data->GetChi2();
    m_chi2 += K40Data->GetChi2();
    m_chi2 += nHData->GetChi2();
    m_chi2 += Co60Data->GetChi2();
    ////m_chi2 += Tl208Data->GetChi2();
    m_chi2 += nC12Data->GetChi2();
    m_chi2 += O16Data->GetChi2();
    m_chi2 += nFe56Data->GetChi2();

    return m_chi2;
}


void gammaResChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


void gammaResChiFunction::SetParameters(double* par)
{
    electronResol::setpA(par[0]);
    electronResol::setpB(par[1]);
    electronResol::setpC(par[2]);
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
    gammaResMinuit->mnparm(iPar, "pA", 0.028, 0.0001, 0.020, 0.030, ierrflag); iPar++;
    gammaResMinuit->mnparm(iPar, "pB", 6.5e-3, 1e-6, 6.0e-3, 7.0e-3, ierrflag); iPar++;
    gammaResMinuit->mnparm(iPar, "pC", 0.0, 0.001, 0.0, 0.00, ierrflag); iPar++;
    
    //gammaResMinuit->FixParameter(0);
    //gammaResMinuit->FixParameter(1);
    //gammaResMinuit->FixParameter(2);
    //gammaResMinuit->FixParameter(3);

    // Minimization strategy
    gammaResMinuit->SetErrorDef(1);
    arglist[0]=2;
    gammaResMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 5000; //maxCalls
    arglist[1] = 0.01; // tolerance
    gammaResMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

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

bool gammaResChiFunction::GridSearch()
{
    double m_like = 1000;
    double tmp_init_pA = 0.0260;
    double tmp_init_pB = 0.0068;
    double tmp_init_pC = 0.;

    int tag_pA = 4;
    int tag_pB = 4;
    int tag_pC = 4;

    // parameters limits:
    double pAHigh = 0.040;
    double pALow  = 0.010;
    double pBHigh = 9.0e-3;
    double pBLow  = 5.0e-3;
    double pCHigh = 0.01;
    double pCLow  = 0.00; 

    double delta = 1e10;

    for(int init_val = 0;init_val < 1;init_val++){
        double step_pA = 0.001;
        double step_pB = 1e-4;
        double step_pC = 0.001;
        for(int iteration = 0; iteration < 200; iteration++){
            for(int bin_pA = -1; bin_pA<2; bin_pA++){
                for(int bin_pB=1; bin_pB<2; bin_pB++) {
                    for(int bin_pC=1; bin_pC<2; bin_pC++) {
                        double tmp_pA = tmp_init_pA + bin_pA*step_pA;
                        double tmp_pB = tmp_init_pB + bin_pB*step_pB;
                        double tmp_pC = tmp_init_pC + bin_pC*step_pC;
                        if(tmp_pA<=pAHigh and tmp_pA>=pALow and tmp_pB<=pBHigh and tmp_pB >= pBLow and tmp_pC<=pCHigh and tmp_pC>=pCLow) {
                            double d = delta;
                            double par[3] = {tmp_pA, tmp_pB, tmp_pC};
                            SetParameters(par);
                            d = GetChi2();

                            if(d<delta){
                                tag_pA = bin_pA;
                                tag_pB = bin_pB;
                                tag_pC = bin_pC;
                                delta = d;
                            }
                        } else continue;
                    }
                }
            }
            tmp_init_pA = tmp_init_pA + static_cast<double>(tag_pA*step_pA);
            tmp_init_pB = tmp_init_pB + static_cast<double>(tag_pB*step_pB);
            tmp_init_pC = tmp_init_pC + static_cast<double>(tag_pC*step_pC);
            if(tag_pA==0 and tag_pB==0 and tag_pC==0) {
                step_pA /= 2.;
                step_pB /= 2.;
                step_pC /= 2.;
            }
            tag_pA = tag_pB = tag_pC = 0;
            cout << tmp_init_pA << " " << tmp_init_pB << " " << tmp_init_pC << " " << delta << endl;
            cout << "iteration: " << iteration << " " << delta << endl; 
        }
    }

    if (delta<m_like) {
        final_pA = tmp_init_pA;
        final_pB = tmp_init_pB;
        final_pC = tmp_init_pC;
    }
    m_DoFit = true;
    m_gridSearch = true;
    return true;

}



void gammaResChiFunction::Plot()
{

    if(!m_DoFit) {
        cout << " >>> Fitting has not been finished !<<< " << endl; return;
    }
    cout << " >>> Plot Resolution Fitting Results <<< " << endl;

    if(m_gridSearch) {
        int nPoints = 3;
        double par[nPoints];
        par[0] = final_pA; par[1] = final_pB; par[2] = final_pC;
        SetParameters(par);
    } else {
        int nPoints = m_nParameter;
        double par[nPoints];
        for(int iPoint=0; iPoint<nPoints; iPoint++) {
            par[iPoint] = m_bestFit[iPoint];
        }
        SetParameters(par);
    }

    TGraphErrors* gData = new TGraphErrors();
    TGraph* gCalc = new TGraph();

    int index = 0;
    for(int iData=0; iData<m_nData; iData++) {
        std::string source = source_name[iData];
        gammaResol* tmpGammaData = mapGammaResol.find(source)->second;
        tmpGammaData->calcGammaNPE();
        double tmp_E       = tmpGammaData->GetEvis();
        double tmp_pred    = tmpGammaData->GetResolPred();
        double tmp_data    = tmpGammaData->GetResolData();
        double tmp_dataErr = tmpGammaData->GetResolDataErr(); 
        cout << tmp_E << " " << tmp_data << " " << tmp_pred << endl;
        gData->SetPoint(index, tmp_E, tmp_data);
        gData->SetPointError(index, 0, tmp_dataErr);
        gCalc->SetPoint(index, tmp_E, tmp_pred);
        index++;
    }

    gData->SetMarkerStyle(20);
    gData->SetMarkerColor(kBlue+1);
    gData->SetLineColor(kBlue+1);
    gData->SetLineWidth(3);
    gData->SetMarkerSize(1.2);
    gCalc->SetMarkerStyle(21);
    gCalc->SetMarkerColor(kPink+2);
    gCalc->SetMarkerSize(1.2);
    gCalc->SetLineColor(kPink+2);
    gCalc->SetLineWidth(3);

    TCanvas* cc = new TCanvas();
    cc->cd(); cc->SetGrid();
    gData->SetTitle("Resolution Fitting; Evis/MeV; Resolution");
    gData->GetYaxis()->SetRangeUser(0.01,0.045);
    gData->Draw("APL");
    gCalc->Draw("PL SAME");
    TLegend* led = new TLegend();
    led->SetFillColor(kWhite);
    led->AddEntry(gData, "data", "PL");
    led->AddEntry(gCalc, "calc", "PL");
    led->Draw("SAME");

    TFile* file = new TFile("./output/gamma/gammaResFit.root", "recreate");
    cc->Write();
    file->Close();
}


