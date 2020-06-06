void read_dist()
{
    TH1D* h1 = new TH1D("h1","", 600, 0, 12);
    TFile* file = TFile::Open("../../fit_param/data/primary_nC12.root");
    TTree* evt = (TTree*)file->Get("evt");
    std::vector<double>* m_elecKE=0; TBranch* b_elecKE;
    evt->SetBranchAddress("ElectronKE", &m_elecKE, &b_elecKE);
    for(int i=0; i<evt->GetEntries(); i++){
        evt->GetEntry(i);
        for(int j=0; j<(*m_elecKE).size(); j++) {
            h1->Fill((*m_elecKE)[j]);
        }
    }
    for(int i=0; i<h1->GetNbinsX();i++) {
        cout << h1->GetBinCenter(i+1) << " " << h1->GetBinContent(i+1) << endl;
    }
}
