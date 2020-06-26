// this script is for pre-sampling for energy resolution model

#include <iostream>
#include <fstream>

using namespace std;

double elecEtrue[810];
double elecNonl[810];
double elecTotPE[810];
double elecTotPESigma[810];
double etrue[9] = {0.662, 0.834, 1.461, 2.223, 2.505, 2.614, 4.945, 6.13, 7.637};
double scale = 3235/2.226;

void load_nonl(string file);
void calc_nonl();
void load_resol(string file);



void pre_sample()
{
    string nonl = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/MC_NonL.txt";
    string resol = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/totPE_smearing.txt";
    string output = "../data/naked_gamma/SecElecNumber.root";

    string input[1];
    for(int i=0; i<1; i++) {
        //input[0] = "Cs137";
        //input[1] = "Mn54";
        //input[0] = "K40";
        //input[0] = "nH";
        //input[0] = "Co60";
        //input[0] = "Tl208";
        //input[0] = "nC12";
        input[0] = "O16";
        //input[0] = "nFe56";
    }

    calc_nonl();
    load_resol(resol);

    Double_t mean[5000]; Double_t sigma[5000]; Double_t SecTotE[5000];  // 1000 gaussian distributions
    TH1D* h1 = new TH1D("Cs137","", 4000, 5000, 12000);  // energy distributions
    TH1D* h2 = new TH1D("Cs137Num", "", 50, 0, 100);  // secondary particle numbers
    TGraph* g1 = new TGraph();

    vector<double> EprmElec;
    ifstream in;
    //TFile* out = TFile::Open(output.c_str(), "RECREATE");

    for(int iFile=0; iFile<1; iFile++) {

        for(int i=0; i<5000; i++) {
            mean[i] = 0; sigma[i] = 0; SecTotE[i] = 0;
        }

        h1->Reset();
        h2->Reset();
        h1->SetName(input[iFile].c_str());
        h2->SetName((input[iFile]+"Num").c_str());

        string filename = "../data/naked_gamma/"+input[iFile]+".txt";
        in.open(filename.c_str());
        if(!in) {cout << "No such input file !!" << endl; break;}
        string line;
        Double_t tmp_elec; Int_t index = 0; int num = 0;  double mean_num = 0;
        while(getline(in,line)) {
            //cout << "Processing " << index << " event" << endl;
            num = 0;
            istringstream ss(line);
            while(ss>>tmp_elec) {
                EprmElec.push_back(tmp_elec); num++;
            }
            mean_num += num;
            h2->Fill(num);
            //double tmp_Evis = 0; double tmp_sigma = 0;
            for(int j=0; j<num; j++) {
                int idx;
                if(EprmElec[j]<0.01) { idx = int(EprmElec[j]/0.001)+1; }
                else {idx = int(EprmElec[j]/0.01)+10;}
                mean[index] += EprmElec[j] * elecNonl[idx]* scale;
                sigma[index] += elecTotPESigma[idx]*elecTotPESigma[idx];
            }

            sigma[index] = TMath::Sqrt(sigma[index]);
            cout << mean[index] << " " << sigma[index] << endl ;
            index++;
            EprmElec.clear();

            if (index==5000) break;
        }

        //g1->SetPoint(iFile, etrue[iFile], mean_num/1000);

        // 2-layer sampling:
        int times = 100000;
        for(int iTime=0; iTime<times; iTime++) {
            int dist = int(gRandom->Uniform(0,5000));
            double evis_tmp = gRandom->Gaus(mean[dist], sigma[dist]);
            h1->Fill(evis_tmp);
        }

        //h2->Write();
        //h1->Write();
        in.close();

    }

    //out->Close();

    //h1->SetTitle(";all primary e+- deposit energy/MeV; entries");
    h1->Draw();

    //g1->SetMarkerColor(kPink+2);
    //g1->SetLineColor(kPink+2);
    //g1->SetLineWidth(2);
    //g1->SetMarkerStyle(20);
    //g1->SetMarkerSize(1);
    //g1->GetYaxis()->SetRangeUser(15,19);
    //g1->GetYaxis()->SetTitle("primary e+- number");
    //g1->GetXaxis()->SetTitle("true gamma energy/MeV");
    //g1->Draw("APL");
}


void load_nonl(string file) {
    ifstream in; in.open(file.c_str());
    if(!in) {cout << "No nonlinearity file !!" << endl;} 
    string line; double tmp_nonl, tmp_E; int num = 0;
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> tmp_E >> tmp_nonl;
        elecEtrue[num] = tmp_E; elecNonl[num] = tmp_nonl*1481.06/scale; num++;
    }
    in.close();

}


void load_resol(string file) {
    ifstream in; in.open(file.c_str());
    if(!in) {cout << "No resolution file !!" << endl;} 
    string line; int num = 0;
    double tmp_E, tmp_totpe, tmp_sigma;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_E >> tmp_totpe >> tmp_sigma;
        elecTotPE[num] = tmp_totpe; elecTotPESigma[num] = tmp_sigma; num++;
    }
    in.close();

}
void calc_nonl()
{

    double p0 = 1.02561e+00;
    double p1 = 1.12245e-01;
    double p2 = 1.39421e+00;
    double p3 = 5.55117e-04;
    for(int i=0; i<800; i++) {
        double eTrue = i/100.;
        double tmp_nonl = (p0+p3*eTrue)/(1+p1*TMath::Exp(-p2*eTrue));
        elecEtrue[i] = eTrue;
        elecNonl[i] = tmp_nonl;
    }
}
