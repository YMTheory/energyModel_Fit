void validate()
{

    TGraph* gCalc;
    TFile* file = TFile::Open("./predSpec.root");
    gCalc = (TGraph*)file->Get("predSpec");
    gCalc->SetMarkerColor(kBlue);
    gCalc->SetLineColor(kBlue);
    gCalc->SetTitle("B12 Spectrum; betaKE/MeV; a.u.");
    gCalc->Draw("AP");

    TH1D* hSimul = new TH1D("Simul", "", 200, 0, 13.4);
    ifstream in;
    in.open("../../B12_edep.txt");
    string line; double tmp;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmp;
        hSimul->Fill(tmp);
    }
    double scale2 = hSimul->GetMaximum();
    hSimul->Scale(1/scale2*1.05);
    hSimul->SetLineColor(kMagenta);
    hSimul->Draw("SAME");

    TLegend* led = new TLegend();
    led->AddEntry(gCalc, "calculatio", "l");
    led->AddEntry(hSimul, "gendacay simulation", "l");
    led->Draw("SAME");
}
