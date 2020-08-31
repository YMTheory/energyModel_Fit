#include "TMinuit.h"

double scale = 3350/2.22;
double n = 2.7e29;
double T = 56.4e-6;
double R = 10e-9;
double normalize = 0.167;
double mass = 0.511;
double chi = 5.1e-29;
double etrue[20];
double data_NL[20];
double data_NL_Err[20];
double pred_NL[20];

void SetParameters(double *par);
double GetChi2();
void LoadData();
void CalcTheory();
double GetChiSquare();
void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

void quench_theory_fit()
{
    // Load Nonlinearity Data 
    LoadData();

    // Theory Calculation

    GetChiSquare();

    TGraphErrors* gNonlData = new TGraphErrors();
    gNonlData->SetLineColor(kRed+1);
    gNonlData->SetLineWidth(2);
    gNonlData->SetMarkerColor(kRed+1);
    gNonlData->SetMarkerStyle(20);
    gNonlData->SetMarkerSize(1);
    TGraph *gNonlCalc = new TGraph();
    gNonlCalc->SetLineColor(kBlue+1);
    gNonlCalc->SetLineWidth(3);
    gNonlCalc->SetMarkerColor(kBlue+1);
    gNonlCalc->SetMarkerStyle(20);
    gNonlCalc->SetMarkerSize(1);
    TLegend* led = new TLegend();
    led->AddEntry(gNonlData, "data", "PL");
    led->AddEntry(gNonlCalc, "theory", "L");

    for(int i=0; i<20; i++)
    {
        gNonlData->SetPoint(i, etrue[i], data_NL[i]);
        gNonlData->SetPointError(i, 0, data_NL_Err[i]);
        gNonlCalc->SetPoint(i, etrue[i], pred_NL[i]);
    }

    gNonlData->SetTitle("Electron Quenching Nonlinearity; etrue/MeV; Nonlinearity");
    gNonlData->Draw("AP");
    gNonlCalc->Draw("L SAME");
    led->Draw("SAME");

}

void SetParameters(double* par) 
{
    n = par[0];
    T = par[1];
    R = par[2];
    normalize = par[3];
}

void LoadData()
{
    ifstream in;
    in.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/totPE_smearing.txt");
    string line; Double_t tmp_etrue = 0; Double_t tmp_pe = 0; Int_t index = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmp_etrue >> tmp_pe;
        if(tmp_etrue < 20) continue;
        tmp_etrue /= 1000;
        etrue[index] = tmp_etrue;
        data_NL[index] = tmp_pe/scale/tmp_etrue;
        data_NL_Err[index] = tmp_pe/scale/tmp_etrue*0.005;
        index++; if(index==20) break;
    }
    
}

void CalcTheory()
{

    double binResol = 0.0001;
    for(int iPoint=0; iPoint<20; iPoint++) {
        double sum = 0;
        for(double inteE=binResol; inteE+binResol<etrue[iPoint]; inteE+=binResol) {
            double gamma_low = (inteE+mass)/mass;
            double beta2_low = 1 - 1/gamma_low/gamma_low;
            double dEdx_low = chi/beta2_low*n*(TMath::Log(inteE*inteE*(1+gamma_low)/2/T/T)-beta2_low*beta2_low);
            double sigma_low = chi/beta2_low/T*(TMath::Log(2*mass*mass*beta2_low*gamma_low*gamma_low/T)-beta2_low);
            double exp_low = TMath::Exp(-n*R*sigma_low);
            double value_low = sigma_low/(dEdx_low)*exp_low;

            double gamma_high = (inteE+binResol+mass)/mass;
            double beta2_high = 1 - 1/gamma_high/gamma_high;
            double dEdx_high = chi/beta2_high*n*(TMath::Log((inteE+binResol)*(inteE+binResol)*(1+gamma_high)/2/T/T)-beta2_high*beta2_high);
            double sigma_high = chi/beta2_high/T*(TMath::Log(2*mass*mass*beta2_high*gamma_high*gamma_high/T)-beta2_high);
            double exp_high = TMath::Exp(-n*R*sigma_high);
            double value_high = sigma_high/(dEdx_high)*exp_high;

            double area = (value_low+value_high)*binResol/2.;
            sum += area;
        }
    
        pred_NL[iPoint] = sum/etrue[iPoint]*normalize*n/scale;
        //gNonlCalc->SetPoint(iPoint, etrue[iPoint], pred_NL[iPoint]);
    }
    
}

double GetChi2()
{
    CalcTheory();
    double m_chi2 = 0;
    for(int i=0; i<20; i++) {
        m_chi2 += TMath::Power((pred_NL[i]-data_NL[i])/data_NL_Err[i], 2);
    }
    return m_chi2;
}

double GetChiSquare()
{
    TMinuit* minuit = new TMinuit();
    minuit->SetFCN(ChisqFCN);
    minuit->SetPrintLevel(1);

    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    minuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    minuit->mnparm(iPar, "n", 2.3e29, 1e28, 0, 0, ierrflag); iPar++;
    minuit->mnparm(iPar, "T", 56.4e-6, 1e-7, 0, 0, ierrflag); iPar++;
    minuit->mnparm(iPar, "R", 10e-9, 1e-9, 0, 0, ierrflag);iPar++;
    minuit->mnparm(iPar, "norm", 0.167, 0.0001, 0, 0, ierrflag); iPar++;


    minuit->SetErrorDef(1);
    arglist[0] = 2;
    minuit->mnexcm("SET STR", arglist, 1, ierrflag);
     
    arglist[0] = 5000;
    arglist[1] = 1;
    minuit->mnexcm("MIGRAD", arglist, 2, ierrflag);

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    minuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;
    delete minuit;
    return min;


}



void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}













