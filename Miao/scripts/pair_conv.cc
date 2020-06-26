void pair_conv()
{
    gStyle->SetOptStat(0);

    TH1D* hSecTotE_conv = new TH1D("hSecTotE_conv", "", 100, 7.625, 7.65);

    string filename = "../data/naked_gamma/nFe56_conv.txt";
    ifstream in;
    in.open(filename.c_str());

    std::vector<double> EprmElec;
    Double_t secTotE[4000];
    Int_t index = 0; string line;
    Double_t tmp_elec;
    for(int i=0; i<4000; i++) {
        secTotE[i] = 0;
    }

    while(getline(in,line)) {
        //cout << "Processing " << index << " event" << endl;
        int num = 0;
        istringstream ss(line);
        while(ss>>tmp_elec) {
            EprmElec.push_back(tmp_elec); num++;
        }
        //double tmp_Evis = 0; double tmp_sigma = 0;
        for(int j=0; j<num; j++) {
            //int idx;
            //if(EprmElec[j]<0.01) { idx = int(EprmElec[j]/0.001)+1; }
            //else {idx = int(EprmElec[j]/0.01)+10;}
            secTotE[index] = secTotE[index]+EprmElec[j];
        }
        hSecTotE_conv->Fill(secTotE[index]); 
        index++;
        EprmElec.clear();

        if (index==4000) break;
    }

    hSecTotE_conv->Draw();
    in.close();



    // other process ...
    
    TH1D* hSecTotE = new TH1D("hSecTotE", "",100, 7.625, 7.65);

    filename = "../data/naked_gamma/nFe56_noconv.txt";
    in.open(filename.c_str());

    for(int i=0; i<4000; i++) {
        secTotE[i] = 0;
    }
    index = 0;

    while(getline(in,line)) {
        //cout << "Processing " << index << " event" << endl;
        int num = 0;
        istringstream ss(line);
        while(ss>>tmp_elec) {
            EprmElec.push_back(tmp_elec); num++;
        }
        //double tmp_Evis = 0; double tmp_sigma = 0;
        for(int j=0; j<num; j++) {
            //int idx;
            //if(EprmElec[j]<0.01) { idx = int(EprmElec[j]/0.001)+1; }
            //else {idx = int(EprmElec[j]/0.01)+10;}
            secTotE[index] = secTotE[index]+EprmElec[j];
        }
        hSecTotE->Fill(secTotE[index]); 
        index++;
        EprmElec.clear();

        if (index==4000) break;
    }

    hSecTotE->GetXaxis()->SetTitle("sum of primary e+- energy/ MeV");
    //hSecTotE->GetXaxis()->SetLimits(0, 2.507);
    hSecTotE->SetLineWidth(3);
    hSecTotE->SetLineColor(kPink+1);
    hSecTotE_conv->SetLineColor(kBlue);
    hSecTotE_conv->SetLineWidth(3);
    hSecTotE->Draw();
    hSecTotE_conv->Draw("SAME");
    TLegend* led = new TLegend();
    led->AddEntry(hSecTotE, "other processes", "l");
    led->AddEntry(hSecTotE_conv, "pair converesion", "l");
    led->Draw("SAME");
    TGraph* g1 = new TGraph();
    g1->SetPoint(0, 7.637, 0);
    g1->SetPoint(1, 7.637, 500);
    g1->SetLineColor(kRed+1);
    g1->SetLineWidth(4);
    g1->Draw("L SAME");
    TText* tex = new TText(6.121, 300, "Etrue for O16");
    tex->Draw("SAME");

}
