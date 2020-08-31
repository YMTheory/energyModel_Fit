
#include <iostream>
#include <fstream>

using namespace std;

double elecEtrue[810];
double elecNonl[810];
double elecTotPE[810];
double elecTotPESigma[810];
double etrue[9] = {0.662, 0.834, 1.461, 2.223, 2.505, 2.614, 4.945, 6.13, 7.637};
//double evis[9]  = {0.606, 0.775, 1.4043, 2.18, 2.470, 2.583, 4.984, 6.208, 7.765};
double evis[9] = {0.61777521, 0.78993756, 1.4312457, 2.2218714, 2.5170658, 2.6319784, 5.0795462, 6.3270696, 7.9131376};
double evisErr[9] = {0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.008, 0.008, 0.014};
double restrue[9] = {0.036244152, 0.032761324, 0.026379808, 0.022626200, 0.021339530, 0.021111111, 0.015659713, 0.014127243, 0.012495652};
double resErr[9] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
double totPE[9] = {897.8, 1148, 2080, 3229, 3658, 3825, 7382, 9195, 11500};
double totPESigma[9] = {32.54, 37.61, 54.87, 73.06, 78.06, 80.75, 115.6, 129.9, 143.7};
double scale = 3350/2.220;

void load_nonl(string file);
void calc_nonl();
void load_resol(string file);
double interpolate_nonl(int idx, double E);
double interpolate_resol(int idx, double E);



void gamma_check()
{
    
    gRandom->SetSeed(44);
    
    string nonl = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/MC_NonL.txt";
    string resol = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/totPE_smearing.txt";
    string output = "../data/naked_gamma/SecElecNumber.root";

    const int nFiles = 9;
    string input[nFiles];
    for(int i=0; i<nFiles; i++) {
        input[0] = "Cs137";
        input[1] = "Mn54";
        input[2] = "K40";
        input[3] = "nH";
        input[4] = "Co60";
        input[5] = "Tl208";
        input[6] = "nC12";
        input[7] = "O16";
        input[8] = "nFe56";
    }

    //calc_nonl();
    load_nonl(nonl);
    load_resol(resol);

    Double_t secTotE[5000]; Int_t secTot[5000];
    Double_t mean[5000]; Double_t sigma[5000];   // 1000 gaussian distributions
    TH1D* h1 = new TH1D("Cs137","", 100, 7.6,7.7);  // energy distributions
    TH1D* h2 = new TH1D("Cs137Num", "", 50, 0, 100);  // secondary particle numbers
    TGraph* g1 = new TGraph(); g1->SetName("pred");
    TGraphErrors* g2 = new TGraphErrors(); g2->SetName("data");
    //TGraph* g2 = new TGraph();
    TGraphErrors* g3 = new TGraphErrors();
    TGraph* g4 = new TGraph(); g4->SetName("pred");
    TGraphErrors* g5 = new TGraphErrors(); g5->SetName("data");
    TGraphErrors* g6 = new TGraphErrors();
    TProfile* gProf = new TProfile("sample","",50, 0, 10000 , 0.52, 0.72);
    TGraph* g7 = new TGraph();
    TGraph* g10 = new TGraph();
    TH1D* h3 = new TH1D("h3", "", 20000, 0, 13000);
    TH2D* h4 = new TH2D("h4", "", 50,0,50,100, 4.935, 4.945);
    TGraph *g8 = new TGraph();
    TGraph* g9 = new TGraph();
    TGraph* gOne = new TGraph();
    TGraph* g11 = new TGraph();
    TGraph* g12 = new TGraph();

    vector<double> EprmElec;
    ifstream in;
    //TFile* out = TFile::Open(output.c_str(), "RECREATE");

    for(int iFile=0; iFile<nFiles; iFile++) {

        for(int i=0; i<5000; i++) {
            secTotE[i] = 0.; secTot[i] = 0;
            mean[i] = 0; sigma[i] = 0;
        }

        h1->Reset();
        h2->Reset();
        h3->Reset();
        h1->SetName(input[iFile].c_str());
        h2->SetName((input[iFile]+"Num").c_str());

        string filename = "../data/naked_gamma/"+input[iFile]+"_all.txt";
        in.open(filename.c_str());
        if(!in) {cout << "No such input file !!  "<< filename  << endl;}
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
            secTot[index] = num;
            h2->Fill(num);
            double tmp_Evis = 0; double tmp_sigma = 0;
            for(int j=0; j<num; j++) {
                int idx; 
                if(EprmElec[j]<0.01) { idx = int(EprmElec[j]/0.001)+1; }
                else {idx = int(EprmElec[j]/0.01)+10;}
                double nonl;
                if(idx==0) {mean[index] += EprmElec[j] * elecNonl[idx] * scale; tmp_Evis += EprmElec[j]*elecNonl[idx]; nonl=elecNonl[idx]; }
                else {mean[index] += EprmElec[j] * interpolate_nonl(idx, EprmElec[j]) * scale; tmp_Evis+= EprmElec[j]*interpolate_nonl(idx, EprmElec[j]); nonl=interpolate_nonl(idx, EprmElec[j]); }
                //mean[index] += EprmElec[j] * elecNonl[idx] * scale;
                //int idx1;
                //idx1 = EprmElec[j]/0.01;
                sigma[index] += interpolate_resol(idx, EprmElec[j]) * interpolate_resol(idx, EprmElec[j]) ;
                //sigma[index] += elecTotPESigma[idx]*elecTotPESigma[idx] ;
                //sigma[index] += EprmElec[j]*EprmElec[j]*elecTotPESigma[idx]*elecTotPESigma[idx]/elecTotPE[idx]/elecTotPE[idx];
            }

            //h1->Fill(secTotE[index]);
            //if(iFile==4) cout <<secTotE[index] << " " << tmp_Evis << " " << mean[index] << endl;
            sigma[index] = TMath::Sqrt(sigma[index]);
            //if(iFile==4) cout << index << " " << mean[index] << " " << sigma[index] << endl;
            //if(iFile==2) g7->SetPoint(index, mean[index], sigma[index]);
            //if(iFile==3) g10->SetPoint(index, mean[index], sigma[index]);
            //if(iFile==2 and secTotE[index]>2.505) std::cout << index << " " << secTotE[index] << endl;
            //if(iFile==8) {  h3->Fill(sigma[index]); h4->Fill(secTot[index], secTotE[index]);}
            index++;
            EprmElec.clear();

            if (index==5000) break;
        }

        // calculate smearing from mean value dist
        double m_mean = 0; double m_std = 0;
        for(int i=0; i<2000; i++) {
            m_mean += mean[i];
        } m_mean/=2000;
        for(int i=0; i<2000; i++) {
            m_std += (mean[i]-m_mean) * (mean[i]-m_mean);
        }
        m_std = TMath::Sqrt(m_std/(1999));
        g10->SetPoint(iFile, etrue[iFile], m_std/m_mean);

        //Double_t mean_secTotE = 0;
        //for(int j=0; j<1000; j++) {
        //    mean_secTotE += secTotE[j];
        //}
        //cout << input[iFile] << "  mean secTotE: " << mean_secTotE/1000./etrue[iFile] << endl;
        //g2->SetPoint(iFile, etrue[iFile], mean_num/1000.);

        // 2-layer sampling:
        const int times = 100000;
        double* sample_E = new double[times];
        Double_t Esample = 0 ; 
        for(int iTime=0; iTime<times; iTime++) {
            int dist = int(gRandom->Uniform(0,index));
            //if(mean[dist]<1400 and iFile==2) continue;
            double evis_tmp = gRandom->Gaus(mean[dist], sigma[dist]);
            Esample += evis_tmp;
            sample_E[iTime] = evis_tmp;
            h3->Fill(evis_tmp); 
            //g1->SetPoint(iTime, iTime, evis_tmp);
            //gProf->Fill(iTime, evis_tmp);
        }
        Esample /= times;
        double tmp_sigma = 0;
        for(int i=0; i<times; i++) {
            tmp_sigma += (sample_E[i]-Esample) * (sample_E[i]-Esample);
        }
        tmp_sigma = TMath::Sqrt(tmp_sigma / (times-1) );
        delete []sample_E;

        h3->Fit("gaus", "Q");
        TF1* func = (TF1*)h3->GetFunction("gaus");
        Esample = func->GetParameter(1);
        tmp_sigma = func->GetParameter(2);
        //cout << input[iFile] << " " << etrue[iFile] << " "<< scale << " " << func->GetParameter(1) << " " << func->GetParameter(2) << " "<< endl;
        //cout << etrue[iFile] << " " << totPE[iFile]/scale/etrue[iFile] << " " << func->GetParameter(1)/scale/etrue[iFile] << endl;
        cout << Esample/scale << "," << tmp_sigma/Esample << endl;
        func->Delete();

        g1->SetPoint(iFile, etrue[iFile], Esample/scale/etrue[iFile]);
        g2->SetPoint(iFile, etrue[iFile],totPE[iFile]/scale/etrue[iFile]);  
        g2->SetPointError(iFile, 0, totPE[iFile]/scale/etrue[iFile]*0.005);
        g3->SetPoint(iFile, etrue[iFile], Esample/totPE[iFile]);
        g3->SetPointError(iFile, 0, Esample/scale/etrue[iFile]/(totPE[iFile]/scale/etrue[iFile]*(1-0.005))-Esample/totPE[iFile] );
        g4->SetPoint(iFile, evis[iFile], tmp_sigma/Esample);
        g5->SetPoint(iFile, evis[iFile], totPESigma[iFile]/totPE[iFile]);
        g5->SetPointError(iFile, 0, totPESigma[iFile]/totPE[iFile]*0.01);
        g6->SetPoint(iFile, evis[iFile], tmp_sigma/totPESigma[iFile]);
        g6->SetPointError(iFile, 0, tmp_sigma/Esample/(totPESigma[iFile]/totPE[iFile]*0.99) - tmp_sigma/Esample/(totPESigma[iFile]/totPE[iFile]));
        g8->SetPoint(iFile, evis[iFile], totPE[iFile]);
        g9->SetPoint(iFile, evis[iFile], Esample);
        
        int idx; 
        if(evis[iFile]<0.01) { idx = int(evis[iFile]/0.001)+1; }
        else {idx = int(evis[iFile]/0.01)+10;}
        g11->SetPoint(iFile, evis[iFile], elecTotPESigma[idx]/elecTotPE[idx]);
        g12->SetPoint(iFile, evis[iFile], TMath::Sqrt( totPESigma[iFile]/totPE[iFile]*totPESigma[iFile]/totPE[iFile] - (elecTotPESigma[idx]/elecTotPE[idx])*(elecTotPESigma[idx]/elecTotPE[idx]) ));
        //cout << evis[iFile] << " " << elecTotPESigma[idx]/elecTotPE[idx] << endl;

        //h2->Write();
        in.close();

    }

    // Plotting

    gOne->SetPoint(0,0,1);
    gOne->SetPoint(1,10,1);
    gOne->SetLineColor(kBlue+1);
    gOne->SetLineWidth(2);
    gOne->SetMarkerColor(kBlue+1);
    gOne->SetMarkerStyle(20);
    gOne->SetMarkerSize(1);

    g1->SetMarkerColor(kOrange+1);
    g1->SetLineColor(kOrange+1);
    g1->SetLineWidth(3);
    g1->SetMarkerStyle(21);
    g1->SetMarkerSize(1.2);
    ////g1->GetYaxis()->SetRangeUser(0.91,1.03);
    g1->SetTitle("Gamma Nonlinearity");
    g1->GetYaxis()->SetTitle("nonlinearity");
    g2->SetMarkerColor(kBlue);
    g2->SetLineColor(kBlue);
    g2->SetLineWidth(2);
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(1.2);
    TLegend* le = new TLegend();
    le->AddEntry(g1, "model", "PL");
    le->AddEntry(g2, "data", "PL");
    g3->SetMarkerColor(kPink+2);
    g3->SetLineColor(kPink+2);
    g3->SetLineWidth(3);
    g3->SetMarkerStyle(20);
    g3->SetMarkerSize(1);
    g3->GetYaxis()->SetRangeUser(0.98,1.020);
    g3->GetYaxis()->SetTitle("pred / true");
    g3->GetXaxis()->SetTitle("gamma total edep/MeV");
    g3->GetXaxis()->SetLabelFont(43);
    g3->GetXaxis()->SetLabelSize(18);
    g3->GetXaxis()->SetTitleFont(43);
    g3->GetXaxis()->SetTitleSize(18);
    g3->GetXaxis()->SetTitleOffset(-2.0);
    g3->GetYaxis()->SetLabelFont(43);
    g3->GetYaxis()->SetLabelSize(8);
    g3->GetYaxis()->SetTitleFont(43);
    g3->GetYaxis()->SetTitleSize(15);
    g3->GetYaxis()->SetTitleOffset(1.5);
    TGraph* zone1 = new TGraph();
    zone1->SetPoint(0,0,1+0.0040);
    zone1->SetPoint(1,8,1+0.0040);
    zone1->SetPoint(2,8,1-0.0040);
    zone1->SetPoint(3,0,1-0.0040);
    zone1->SetFillStyle(3005);
    zone1->SetFillColor(kTeal);
    TText *text1 = new TText(3, 1.0042, "0.4% relative bias zone");
    text1->SetTextFont(43);
    text1->SetTextSize(18);

    TCanvas *c1 = new TCanvas("nonl",
       "nonlinearity model for positron",10,40,800,600);
    c1->Range(0,0,25,18);
    c1->SetFillColor(0);
    TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.12,0.98,0.40,0);
    TPad *pad2 = new TPad("pad2","This is pad2",0.02,0.345,0.98,0.93,0);
    pad2->Draw();
    pad1->Draw();

    pad1->cd();
    g3->Draw("AP");
    zone1->Draw("F same");
    text1->Draw("SAME");
    gOne->Draw("L SAME");

    pad2->cd();
    g1->Draw("APL");
    g2->Draw("P SAME");
    le->Draw("SAME");
    

    g4->SetMarkerColor(kOrange+1);
    g4->SetLineColor(kOrange+1);
    g4->SetLineWidth(3);
    g4->SetMarkerStyle(21);
    g4->SetMarkerSize(1.2);
    g5->SetMarkerColor(kBlue);
    g5->SetLineColor(kBlue);
    g5->SetLineWidth(2);
    g5->SetMarkerStyle(20);
    g5->SetMarkerSize(1.2);
    g11->SetMarkerColor(kViolet+1);
    g11->SetLineColor(kViolet+1);
    g11->SetLineWidth(3);
    g11->SetMarkerStyle(21);
    g11->SetMarkerSize(1.2);
    g4->SetTitle("Gamma Resolution");
    g4->GetYaxis()->SetTitle("Resolution");
    TLegend* ld = new TLegend();
    ld->AddEntry(g4, "model", "PL");
    ld->AddEntry(g5, "data", "PL");
    ld->Draw("SAME");
    g6->SetMarkerColor(kPink+2);
    g6->SetLineColor(kPink+2);
    g6->SetLineWidth(3);
    g6->SetMarkerStyle(20);
    g6->SetMarkerSize(1);
    g6->GetYaxis()->SetRangeUser(0.90,1.10);
    g6->GetYaxis()->SetTitle("pred nonl/true nonl");
    g6->GetXaxis()->SetTitle("true gamma energy/MeV");
    g6->GetYaxis()->SetTitle("pred / true");
    g6->GetXaxis()->SetTitle("gamma evis/MeV");
    g6->GetXaxis()->SetLabelFont(46);
    g6->GetXaxis()->SetLabelSize(18);
    g6->GetXaxis()->SetTitleFont(46);
    g6->GetXaxis()->SetTitleSize(18);
    g6->GetXaxis()->SetTitleOffset(-2.0);
    g6->GetYaxis()->SetLabelFont(46);
    g6->GetYaxis()->SetLabelSize(8);
    g6->GetYaxis()->SetTitleFont(46);
    g6->GetYaxis()->SetTitleSize(15);
    g6->GetYaxis()->SetTitleOffset(1.5);
    TGraph* zone2 = new TGraph();
    zone2->SetPoint(0,0,1+0.015);
    zone2->SetPoint(1,8,1+0.015);
    zone2->SetPoint(2,8,1-0.015);
    zone2->SetPoint(3,0,1-0.015);
    zone2->SetFillStyle(3005);
    zone2->SetFillColor(kTeal);
    TText *text2 = new TText(3, 1.016, "1.5% relative bias zone");
    text2->SetTextFont(43);
    text2->SetTextSize(18);

    // Create a new canvas.
    TCanvas *c2 = new TCanvas("resol",
       "resolution model for positron",10,40,800,600);
    c2->Range(0,0,25,18);
    c2->SetFillColor(0);
    TPad *pad3 = new TPad("pad3","This is pad3",0.02,0.12,0.98,0.40,0);
    TPad *pad4 = new TPad("pad4","This is pad4",0.02,0.345,0.98,0.93,0);
    pad4->Draw();
    pad3->Draw();

    pad3->cd();
    g6->Draw("AP");
    zone2->Draw("F same");
    text2->Draw("SAME");
    gOne->Draw("L SAME");

    pad4->cd();
    g4->Draw("APL");
    g5->Draw("P SAME");
    g11->Draw("PL SAME");
    ld->Draw("SAME");



    TCanvas* c5 = new TCanvas();
    g12->SetMarkerColor(kRed+1);
    g12->SetMarkerStyle(20);
    g12->SetMarkerSize(1.3);
    g12->SetLineColor(kRed+1);
    g12->SetLineWidth(2);
    g12->Draw("APL");


    TFile* ff = new TFile("gamma2.root", "recreate");
    g4->Write();
    g5->Write();
    ff->Close();

}


void load_nonl(string file) {
    ifstream in; in.open(file.c_str());
    if(!in) {cout << "No nonlinearity file !!" << endl;} 
    string line; double tmp_E, tmp_nonl; int num = 0;
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
    for(int i=0; i<8000; i++) {
        double eTrue = i/1000.;
        double tmp_nonl = (p0+p3*eTrue)/(1+p1*TMath::Exp(-p2*eTrue));
        elecEtrue[i] = eTrue;
        elecNonl[i] = tmp_nonl;
    }
}

double interpolate_nonl(int idx, double E) {
    double delta = (elecNonl[idx]-elecNonl[idx-1])*(E-elecEtrue[idx-1]/1000)/(elecEtrue[idx]/1000-elecEtrue[idx-1]/1000) ;
    return delta+elecNonl[idx-1];
}


double interpolate_resol(int idx, double E) {
    double delta = (elecTotPESigma[idx]-elecTotPESigma[idx-1])*(E-elecEtrue[idx-1]/1000)/(elecEtrue[idx]/1000-elecEtrue[idx-1]/1000) ;
    return delta+elecTotPESigma[idx-1];
}
