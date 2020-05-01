void test(){
    TH1D* h1; TH1D* h2 = new TH1D("h2", "", 1400, 0, 14);
    TFile* f = TFile::Open("data/electron/B12_file.root");
    h1 = (TH1D*)f->Get("B12Origin");
    //cout << "bin number: " << h1->GetNbinsX() << endl;
    for(int i=0; i<h1->GetNbinsX(); i++) {
        h2->SetBinContent(i, h1->GetBinContent(i));
        cout << i << " " << h1->GetBinContent(i) << endl;
    }
    f->Close();
    cout << "bin number: " << h2->GetNbinsX() << endl;
}
