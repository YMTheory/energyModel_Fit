void test(){
    TGraph* g1 = new TGraph();
    ifstream in;
    in.open("./log");
    string line; double tmp_E, tmp_nl; int i=0;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_E >> tmp_nl;
        g1->SetPoint(i, tmp_E, tmp_nl); i++;
    }
    g1->SetMarkerSize(0.5);
    g1->SetMarkerColor(kBlue+1);
    g1->Draw("APL");

    TFile* file = TFile::Open("./data/Gamma_Electron1.root");
    TH1D* hnFe56 = (TH1D*)file->Get("gammanFe56");
    hnFe56->SetMarkerSize(0.5);
    hnFe56->SetMarkerColor(kRed+1);
    hnFe56->SetLineColor(kRed+1);
    hnFe56->Draw("PEX0 SAME");

}

