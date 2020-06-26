double elecEtrue[810];
double elecNonl[810];
double elecTotPE[810];
double elecTotPESigma[810];
double scale = 3450/2.220;
double gamma_totPE; double gamma_totPESigma;
double gamma_Etrue = 0.511; //MeV
double KEtrue[8] = { 0, 0.1,0.2,0.5,1,2,5,7.98};
double totPE[8] = { 1.36313e+03 ,1.49939e+03,1.64069e+03, 2.08453e+03, 2.83125e+03, 4.38023e+03, 8.99394e+03, 13600};
double totPESigma[8] = { 3.90163e+01, 4.04304e+01, 4.18144e+01, 4.73005e+01, 5.63472e+01,7.25598e+01, 1.13471e+02, 148};

void load_nonl(string file);
void load_resol(string file);
double interpolate_nonl(int idx, double E);
double interpolate_resol(int idx, double E);
void pred_511keVGamma();

void positron_check()
{
    string nonl = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/MC_NonL.txt";
    string resol = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/totPE_smearing.txt";
    load_nonl(nonl);
    load_resol(resol);

    pred_511keVGamma();
    cout << "511keV gamma: " << gamma_totPE/scale << " " << gamma_totPESigma/gamma_totPE << endl;

    double mass = 1.022;  //MeV

    TGraph* gNonlCalc = new TGraph();
    TGraphErrors* gNonlData = new TGraphErrors();
    TGraphErrors* gNonlRatio = new TGraphErrors();
    TGraph* gResCalc = new TGraph();
    TGraphErrors* gResData = new TGraphErrors();
    TGraphErrors* gResRatio = new TGraphErrors();
    TGraph* gOne = new TGraph();

    TH1D* hPositron = new TH1D("hPositron", "", 5000, 0, 15000);

    for(int i=0; i<8; i++) {
        hPositron->Reset();

        // kinetic energy part 
        double ketrue = KEtrue[i];
        double etrue = ketrue+mass;
        int idx;
        if(ketrue<0.01) { idx = int(ketrue/0.001)+1; }
        else {idx = int(ketrue/0.01)+10;}
        double ke_nonl = interpolate_nonl(idx, ketrue);
        double ke_sigma = interpolate_resol(idx, ketrue);  // NPE sigma

        double evis = ketrue * ke_nonl + 2* gamma_totPE/scale;
        double tot_nonl = evis/etrue;

        double tot_sigma = TMath::Sqrt(ke_sigma*ke_sigma + gamma_totPESigma*gamma_totPESigma *2 );
        double tot_resol = tot_sigma / scale / evis;

        for(int i=0; i<100000; i++) {  // 100000 times sampling
            double elec_pe = gRandom->Gaus(elecTotPE[idx], ke_sigma);
            double gamma_pe1 = gRandom->Gaus(gamma_totPE, gamma_totPESigma);
            double gamma_pe2 = gRandom->Gaus(gamma_totPE, gamma_totPESigma);
            hPositron->Fill(elec_pe+gamma_pe1+gamma_pe2);
        }

        hPositron->Fit("gaus", "Q");
        TF1* func = (TF1*)hPositron->GetFunction("gaus");
        double pe_mean = func->GetParameter(1);
        double pe_sigma = func->GetParameter(2);
        func->Delete();

        tot_nonl = pe_mean/scale/etrue;
        tot_resol = pe_sigma / pe_mean;

        gNonlData->SetPoint(i, etrue, totPE[i]/scale/etrue);
        gNonlData->SetPointError(i, 0, totPE[i]/scale/etrue*0.005);
        gNonlCalc->SetPoint(i, etrue, tot_nonl);
        gNonlRatio->SetPoint(i, etrue, tot_nonl/(totPE[i]/scale/etrue));
        gNonlRatio->SetPointError(i, 0, (tot_nonl/(totPE[i]/scale/etrue*0.995)) - (tot_nonl/(totPE[i]/scale/etrue))  );
        gResData->SetPoint(i, totPE[i]/scale, totPESigma[i]/totPE[i]);
        gResData->SetPointError(i, 0, totPESigma[i]/totPE[i]*0.01);
        gResCalc->SetPoint(i, evis, tot_resol);
        gResRatio->SetPoint(i, evis, tot_resol/(totPESigma[i]/totPE[i]));
        gResRatio->SetPointError(i, 0, (tot_resol/(totPESigma[i]/totPE[i]*0.99))-(tot_resol/(totPESigma[i]/totPE[i])) );
    }

    // Plotting
    gOne->SetPoint(0,0,1);
    gOne->SetPoint(1,10,1);
    gOne->SetLineColor(kBlue+1);
    gOne->SetLineWidth(2);
    gOne->SetMarkerColor(kBlue+1);
    gOne->SetMarkerStyle(20);
    gOne->SetMarkerSize(1);

    gNonlCalc->SetLineColor(kBlue+1);
    gNonlCalc->SetLineWidth(2);
    gNonlCalc->SetMarkerColor(kBlue+1);
    gNonlCalc->SetMarkerStyle(20);
    gNonlCalc->SetMarkerSize(1);
    gNonlData->GetYaxis()->SetTitle("Nonlinearity");
    gNonlData->SetTitle("Positron Nonlinearity");
    //gNonlCalc->GetYaxis()->SetRangeUser(0.7, 1.1);
    gNonlData->SetLineColor(kGreen+1);
    gNonlData->SetLineWidth(2);
    gNonlData->SetMarkerColor(kGreen+1);
    gNonlData->SetMarkerStyle(21);
    gNonlData->SetMarkerSize(1);
    gNonlRatio->SetLineColor(kPink+2);
    gNonlRatio->SetLineWidth(2);
    gNonlRatio->SetMarkerColor(kPink+2);
    gNonlRatio->SetMarkerStyle(21);
    gNonlRatio->SetMarkerSize(1);
    gNonlRatio->GetYaxis()->SetRangeUser(0.985, 1.015);
    gNonlRatio->GetYaxis()->SetTitle("pred / true");
    gNonlRatio->GetXaxis()->SetTitle("positron total edep/MeV");
    gNonlRatio->GetXaxis()->SetLabelFont(43);
    gNonlRatio->GetXaxis()->SetLabelSize(18);
    gNonlRatio->GetXaxis()->SetTitleFont(43);
    gNonlRatio->GetXaxis()->SetTitleSize(18);
    gNonlRatio->GetXaxis()->SetTitleOffset(-2.0);
    gNonlRatio->GetYaxis()->SetLabelFont(43);
    gNonlRatio->GetYaxis()->SetLabelSize(8);
    gNonlRatio->GetYaxis()->SetTitleFont(43);
    gNonlRatio->GetYaxis()->SetTitleSize(15);
    gNonlRatio->GetYaxis()->SetTitleOffset(1.5);
    TGraph* zone1 = new TGraph();
    zone1->SetPoint(0,0,1+0.0100);
    zone1->SetPoint(1,9.1,1+0.0100);
    zone1->SetPoint(2,9.1,1-0.0100);
    zone1->SetPoint(3,0,1-0.0100);
    zone1->SetFillStyle(3005);
    zone1->SetFillColor(kTeal);
    TText *text1 = new TText(3, 1.0105, "1.0% relative bias zone");
    text1->SetTextFont(43);
    text1->SetTextSize(18);
    TLegend* l1 = new TLegend();
    l1->AddEntry(gNonlCalc, "calc", "l");
    l1->AddEntry(gNonlData, "data", "l");

    // Create a new canvas.
    TCanvas *c1 = new TCanvas("nonl",
       "nonlinearity model for positron",10,40,800,600);
    c1->Range(0,0,25,18);
    c1->SetFillColor(0);
    TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.12,0.98,0.40,0);
    TPad *pad2 = new TPad("pad2","This is pad2",0.02,0.345,0.98,0.93,0);
    pad2->Draw();
    pad1->Draw();

    pad1->cd();
    gNonlRatio->Draw("AP");
    zone1->Draw("F same");
    text1->Draw("SAME");
    gOne->Draw("L SAME");

    pad2->cd();
    gNonlData->Draw("APL");
    gNonlCalc->Draw("P SAME");
    l1->Draw("SAME");
    

    gResCalc->SetLineColor(kBlue+1);
    gResCalc->SetLineWidth(2);
    gResCalc->SetMarkerColor(kBlue+1);
    gResCalc->SetMarkerStyle(20);
    gResCalc->SetMarkerSize(1);
    gResData->GetYaxis()->SetTitle("Resolution");
    gResData->SetTitle("Positron Resolution");
    //gResCalc->GetYaxis()->SetRangeUser(0.7, 1.1);
    gResData->SetLineColor(kGreen+1);
    gResData->SetLineWidth(2);
    gResData->SetMarkerColor(kGreen+1);
    gResData->SetMarkerStyle(21);
    gResData->SetMarkerSize(1);
    gResRatio->SetLineColor(kPink+2);
    gResRatio->SetLineWidth(2);
    gResRatio->SetMarkerColor(kPink+2);
    gResRatio->SetMarkerStyle(21);
    gResRatio->SetMarkerSize(1);
    gResRatio->GetYaxis()->SetRangeUser(0.90, 1.10);
    gResRatio->GetYaxis()->SetTitle("pred / true");
    gResRatio->GetXaxis()->SetTitle("positron evis/MeV");
    gResRatio->GetXaxis()->SetLabelFont(43);
    gResRatio->GetXaxis()->SetLabelSize(18);
    gResRatio->GetXaxis()->SetTitleFont(43);
    gResRatio->GetXaxis()->SetTitleSize(18);
    gResRatio->GetXaxis()->SetTitleOffset(-2.0);
    gResRatio->GetYaxis()->SetLabelFont(43);
    gResRatio->GetYaxis()->SetLabelSize(8);
    gResRatio->GetYaxis()->SetTitleFont(43);
    gResRatio->GetYaxis()->SetTitleSize(15);
    gResRatio->GetYaxis()->SetTitleOffset(1.5);
    TGraph* zone2 = new TGraph();
    zone2->SetPoint(0,0,1+0.050);
    zone2->SetPoint(1,9.5,1+0.050);
    zone2->SetPoint(2,9.5,1-0.050);
    zone2->SetPoint(3,0,1-0.050);
    zone2->SetFillStyle(3005);
    zone2->SetFillColor(kTeal);
    TText *text2 = new TText(3, 1.052, "5.0% relative bias zone");
    text2->SetTextFont(43);
    text2->SetTextSize(18);
    TLegend* l2 = new TLegend();
    l2->AddEntry(gResCalc, "calc", "l");
    l2->AddEntry(gResData, "data", "l");

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
    gResRatio->Draw("AP");
    zone2->Draw("F same");
    text2->Draw("SAME");
    gOne->Draw("L SAME");

    pad4->cd();
    gResData->Draw("APL");
    gResCalc->Draw("P SAME");
    l2->Draw("SAME");

}


// sub-function definition

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

double interpolate_nonl(int idx, double E) {
    double delta = (elecNonl[idx]-elecNonl[idx-1])*(E-elecEtrue[idx-1]/1000)/(elecEtrue[idx]/1000-elecEtrue[idx-1]/1000) ;
    return delta+elecNonl[idx-1];
}


double interpolate_resol(int idx, double E) {
    double delta = (elecTotPESigma[idx]-elecTotPESigma[idx-1])*(E-elecEtrue[idx-1]/1000)/(elecEtrue[idx]/1000-elecEtrue[idx-1]/1000) ;
    return delta+elecTotPESigma[idx-1];
}

void pred_511keVGamma()
{
    vector<double> EprmElec;
    TH1D* h1 = new TH1D("h511keV", "", 200, 400, 1200);
    ifstream in;
    in.open("../data/naked_gamma/gamma511keV.txt");
    if(!in) {cout << "No such input file !!  "  << endl;}
    string line;
    Double_t tmp_elec; Int_t index = 0; int num = 0;
    double mean[5000]; double sigma[5000];
    for(int i=0; i<5000; i++) {
        mean[i] = 0.; sigma[i] = 0.;
    }

    while(getline(in,line)) {
       num = 0;
       istringstream ss(line);
       while(ss>>tmp_elec) {
           EprmElec.push_back(tmp_elec); num++; 
       }
       double secTotE = 0;
       for(int j=0; j<num; j++) {
           int idx;
           if(EprmElec[j]<0.01) { idx = int(EprmElec[j]/0.001)+1; }
           else {idx = int(EprmElec[j]/0.01)+10;}
           secTotE += EprmElec[j]; 
           if(idx==0) {
           mean[index] += EprmElec[j] * elecNonl[idx] * scale; 
           sigma[index] += elecTotPESigma[idx]*elecTotPESigma[idx]; }
           else {mean[index] += EprmElec[j] * interpolate_nonl(idx, EprmElec[j]) * scale;
                 sigma[index] += interpolate_resol(idx, EprmElec[j]) * interpolate_resol(idx, EprmElec[j]) ;}
        }
        sigma[index] = TMath::Sqrt(sigma[index]);
        index++;
        EprmElec.clear();
        if (index==5000) break;
    }
    // 2-layer sampling:
    const int times = 100000;
    for(int iTime=0; iTime<times; iTime++) {
        int dist = int(gRandom->Uniform(0,index));
        double evis_tmp = gRandom->Gaus(mean[dist], sigma[dist]);
        h1->Fill(evis_tmp);  ;
    }

    h1->Fit("gaus", "Q");
    TF1* func = (TF1*)h1->GetFunction("gaus");
    gamma_totPE = func->GetParameter(1);
    gamma_totPESigma = func->GetParameter(2);
    func->Delete();

    in.close();
}

