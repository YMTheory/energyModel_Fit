double elecEtrue[810];
double elecNonl[810];
double elecTotPE[810];
double elecTotPESigma[810];
double scale = 3350/2.220;
void load_nonl(string file);
void load_resol(string file);

void read()
{
    string nonl = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/MC_NonL.txt";
    string resol = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/totPE_smearing.txt";
    load_nonl(nonl);
    load_resol(resol);

    const int num=6;
    double etru[num] = {0.5, 1, 2, 3, 6, 7};
    TGraph* gElecNonl = new TGraph();
    for(int i=0; i<num; i++) {
        double etrue = etru[i];
        int idx1;
        if(etrue<0.01) { idx1 = int(etrue/0.001)+1; }
        else {idx1 = int(etrue/0.01)+10;}
        double evis = elecTotPE[idx1]/scale;
        double res = elecTotPESigma[idx1]/elecTotPE[idx1];
        //gElecNonl->SetPoint(i, etrue, evis/etrue);
        gElecNonl->SetPoint(i, elecTotPE[idx1]/scale, res);
    }

    TCanvas* cc = new TCanvas();  cc->SetGrid();
    gElecNonl->SetLineColor(kBlue+1);
    gElecNonl->GetXaxis()->SetLimits(0,9.2);
    gElecNonl->GetYaxis()->SetRangeUser(0., 0.040);
    gElecNonl->SetLineWidth(2);
    gElecNonl->SetTitle("; Evis/MeV; resolution");
    gElecNonl->Draw("AL");

    TFile* f1 = TFile::Open("gamma2.root");
    TGraph* gGamCalc = (TGraph*)f1->Get("pred");
    TGraphErrors* gGamData = (TGraphErrors*)f1->Get("data");
    gGamData->SetMarkerColor(kPink+1);
    gGamData->SetMarkerSize(1.0);
    gGamCalc->SetLineColor(kViolet+1);
    gGamCalc->SetLineWidth(2);
    gGamCalc->Draw("L SAME");
    gGamData->Draw("P SAME");

    TFile* f2 = TFile::Open("positron2.root");
    TGraph* gPosiCalc = (TGraph*)f2->Get("pred");
    TGraphErrors* gPosiData = (TGraphErrors*)f2->Get("data");
    gPosiData->SetMarkerColor(kGreen+1);
    gPosiData->SetMarkerSize(1.0);
    gPosiCalc->SetLineColor(kOrange+1);
    gPosiCalc->SetLineWidth(2);
    gPosiCalc->Draw("L SAME");
    gPosiData->Draw("P SAME");

    TLegend* le = new TLegend();
    le->AddEntry(gElecNonl, "electron", "l");
    le->AddEntry(gGamCalc, "gamma", "l");
    le->AddEntry(gPosiCalc, "positron", "l");
    le->Draw("SAME");
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
