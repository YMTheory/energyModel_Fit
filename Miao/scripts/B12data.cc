void B12data()
{
    TFile *file = new TFile("../data/electron/B12.root");
    TH1D* m_edep = (TH1D*)file->Get("Evis");
    m_edep->Draw();

    TFile* newfile = new TFile("B12.root", "recreate");
    m_edep->SetName("b12");
    m_edep->Write();
    newfile->Close();
}
