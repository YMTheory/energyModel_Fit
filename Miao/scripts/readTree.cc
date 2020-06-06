void readTree()
{
    char name[10]; 
    TFile* file = TFile::Open("../data/electron/theo/B12_theo.root");
    TTree* t = (TTree*)file->Get("T");
    t->SetBranchAddress("photonName", name);
    for(int i=0; i<t->GetEntries(); i++) {
        t->GetEntry(i);
        cout <<i << " "  << name[0] << endl;
    }
}
