double elecEtrue[810];
double elecNonl[810];
double elecTotPE[810];
double elecTotPESigma[810];
void load_nonl(string file);
double m_quenchingShape1[1689];
void load_quench(string file);
std::vector<double> m_Etrue;
std::vector<double> m_Cerenkov;
void load_cerenkov(string file);
void load_resol(string file);
double m_energyScale = 3350/2.22;
double etrue[6] = {0.5, 1, 2, 3, 4, 6};
double totPE[6] = {721.551, 1481.06, 3011.85, 4550.84, 6088.3, 9152.69};
double totPESigma[6] = {0.005, 0.005, 0.005, 0.005, 0.005, 0.005};

void electron_check()
{
    string nonl = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/MC_NonL.txt";
    string quenchFile = "/Users/yumiao/Documents/Works/github/energyModel_Fit/Miao/data/electron/Quench5.root";
    string CerenkovFile = "/Users/yumiao/Documents/Works/github/energyModel_Fit/Miao/data/electron/Cer.dat";
    string resol = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/totPE_smearing.txt";

    load_nonl(nonl);
    load_quench(quenchFile);
    load_cerenkov(CerenkovFile);
    load_resol(resol);

    TGraphErrors* g1 = new TGraphErrors();
    TGraph* g2 = new TGraph();
    TGraph* g3 = new TGraph();
    TGraph* g4 = new TGraph();
    TGraph* g5 = new TGraph();

    double kA = 0.962;
    double kC = 1.00;
    double pA = 0.02576;
    double pB = 0.006864;
    double pC = 0;

    //for(int iSample=0; iSample<6; iSample++) {
    //    int idx = 0;
    //    double eTrue = etrue[iSample];
    //    if(eTrue<=0.1) { idx = int(eTrue/0.001); }
    //    else { idx = int((eTrue-0.1)/0.01)+100; }

    //    double quenchNL  =  kA * m_quenchingShape1[idx] ;

    //    int idx1 = 0;
    //    idx1 = eTrue/0.01;
    //    double cerNL     =  kC * m_Cerenkov[idx1];
    //    cerNL = cerNL/m_energyScale/eTrue;

    //    cout << etrue[iSample] << " " << quenchNL << " "<<cerNL << endl;
    //    g1->SetPoint(iSample, etrue[iSample], totPE[iSample]/m_energyScale/etrue[iSample]);
    //    g1->SetPointError(iSample, 0, totPESigma[iSample]*totPE[iSample]/m_energyScale/etrue[iSample]);
    //    g2->SetPoint(iSample, etrue[iSample], quenchNL+cerNL);
    //    g3->SetPoint(iSample, etrue[iSample], (quenchNL+cerNL)/(totPE[iSample]/m_energyScale/etrue[iSample]));
    //
    //}

    for(int i=1; i<80; i++) {
        double eTrue = elecEtrue[i*10]/1000;
        double evis = elecTotPE[i]/m_energyScale;
        int idx = 0;
        if(eTrue<=0.1) { idx = int(eTrue/0.001); }
        else { idx = int((eTrue-0.1)/0.01)+100; }
        double quenchNL  =  kA * m_quenchingShape1[idx] ;

        int idx1 = 0;
        idx1 = eTrue/0.01;
        double cerNL     =  kC * m_Cerenkov[idx1];
        cerNL = cerNL/m_energyScale/eTrue;

        cout << eTrue << " " << elecNonl[i] << " " << quenchNL+cerNL << endl;

        if(eTrue==0.1 or eTrue==0.2 or eTrue==0.3 or eTrue==0.6 or eTrue==0.8) continue;
        g1->SetPoint(i-1, eTrue, elecNonl[i*10]);
        g1->SetPointError(i-1, 0, elecNonl[i*10]*0.005);
        g2->SetPoint(i-1, eTrue, quenchNL+cerNL);
        g3->SetPoint(i-1, eTrue, (quenchNL+cerNL)/(elecNonl[i*10]));
        g4->SetPoint(i-1, evis, elecTotPESigma[i*10]/elecTotPE[i*10]);
        g5->SetPoint(i-1, evis, TMath::Sqrt(pA*pA/evis+pB*pB));

    }

    TCanvas* c1 = new TCanvas(); c1->cd(); c1->SetGrid();
    g1->SetMarkerColor(kOrange+1);
    g1->SetLineColor(kOrange+1);
    g1->SetLineWidth(1);
    g1->SetMarkerStyle(21);
    g1->SetMarkerSize(0.5);
    ////g1->GetYaxis()->SetRangeUser(0.91,1.03);
    g1->SetTitle("electron Nonlinearity");
    g1->GetYaxis()->SetTitle("nonlinearity");
    g1->GetXaxis()->SetTitle("electron etrue/MeV");
    g1->Draw("AP");

    ////TCanvas* c2 = new TCanvas();
    g2->SetMarkerColor(kBlue);
    g2->SetLineColor(kBlue);
    g2->SetLineWidth(2);
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(0.5);
    ////g2->getyaxis()->setrangeuser(25,29);
    ////g2->getyaxis()->settitle("primary e+- number");
    ////g2->getxaxis()->settitle("true gamma energy/mev");
    g2->Draw("PL SAME");
    //g2->Draw("AP");

    TLegend* le = new TLegend();
    le->AddEntry(g1, "data", "PL");
    le->AddEntry(g2, "model", "PL");
    //le->AddEntry(g3, "model average", "PL");
    le->Draw("SAME");

    TCanvas* c2 = new TCanvas();
    g3->SetMarkerColor(kPink+2);
    g3->SetLineColor(kPink+2);
    g3->SetLineWidth(3);
    g3->SetMarkerStyle(21);
    g3->SetMarkerSize(0.8);
    ////g1->GetYaxis()->SetRangeUser(0.91,1.03);
    g3->SetTitle("Electron Nonlinearity");
    g3->GetYaxis()->SetTitle("pred nonl / true nonl");
    g3->GetXaxis()->SetTitle("electron etrue/MeV");
    g3->Draw("AP");
    TGraph* g6 = new TGraph();
    g6->SetPoint(0, 0, 1);
    g6->SetPoint(1, 8, 1);
    g6->SetLineWidth(4);
    g6->SetLineColor(kBlue);
    g6->Draw("L SAME");

    TCanvas* c3 = new TCanvas(); c3->cd(); c3->SetGrid();
    g4->SetMarkerColor(kOrange+1);
    g4->SetLineColor(kOrange+1);
    g4->SetLineWidth(1);
    g4->SetMarkerStyle(21);
    g4->SetMarkerSize(0.5);
    ////g4->GetYaxis()->SetRangeUser(0.91,1.03);
    g4->SetTitle("electron Nonlinearity");
    g4->GetYaxis()->SetTitle("nonlinearity");
    g4->GetXaxis()->SetTitle("electron etrue/MeV");
    g4->Draw("AP");

    ////TCanvas* c2 = new TCanvas();
    g5->SetMarkerColor(kBlue);
    g5->SetLineColor(kBlue);
    g5->SetLineWidth(2);
    g5->SetMarkerStyle(20);
    g5->SetMarkerSize(0.5);
    ////g5->getyaxis()->setrangeuser(25,29);
    ////g5->getyaxis()->settitle("primary e+- number");
    ////g5->getxaxis()->settitle("true gamma energy/mev");
    g5->Draw("P SAME");
    //g5->Draw("AP");

    TLegend* ll = new TLegend();
    ll->AddEntry(g4, "data", "PL");
    ll->AddEntry(g5, "model", "PL");
    //le->AddEntry(g3, "model average", "PL");
    ll->Draw("SAME");


}


void load_nonl(string file) {
    ifstream in; in.open(file.c_str());
    if(!in) {cout << "No nonlinearity file !!" << endl;}
    string line; double tmp_E, tmp_nonl; int num = 0;
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> tmp_E >> tmp_nonl;
        elecEtrue[num] = tmp_E; elecNonl[num] = tmp_nonl*1481.06/m_energyScale; 
        num++;
    }
    in.close();

}

void load_quench(string file)
{
    cout << " >>> Loading Quenching NL Data <<< " << endl;
    TFile* quenchingFile = new TFile(file.c_str(), "read");
    if(!quenchingFile) { std::cout << " >>> Fail to Open QuenchNL File <<< " << std::endl; }

    int kbIdx = 65;
    stringstream ss; ss << kbIdx;
    TString name1 = "kB"+ss.str();
    TH1D* quench1G = (TH1D*)quenchingFile->Get(name1);
    if(!quench1G) { cout << "No Such A Histogram " << name1 << " in Quench.root File" << endl; return;  }
    //double* quench1 = quench1G->GetY();

    for(int sampleIdx=0; sampleIdx<1689; sampleIdx++)
    {
		//m_quenchingShape1[kbIdx][sampleIdx] = quench1[sampleIdx];
        m_quenchingShape1[sampleIdx] = quench1G->GetBinContent(sampleIdx+1);
    }
	quenchingFile->Close();
    delete quenchingFile;

}


void load_cerenkov(string file)
{
    cout << " >>> Loading Electron Cerenkov Shape <<< " << endl;
    ifstream in;
    in.open(file.c_str());
    if(!in){
        cout << " >>> Fail to Open Cerenkov File !! <<< " << endl;
    }
    string line;

    double tmp_Edep, tmp_rel, tmp_abs;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_Edep >> tmp_rel >> tmp_abs ;
        m_Etrue.push_back(tmp_Edep/1000.);
        m_Cerenkov.push_back(tmp_abs);
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
