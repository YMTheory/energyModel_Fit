void convert()
{
    TFile* in = new TFile("B12_theo.root", "read");
    TH1F* hh0 = (TH1F*)in->Get("hh0");
    TH1F* hh1 = (TH1F*)in->Get("hh1");
    TH1F* hh2 = (TH1F*)in->Get("hh2");

    TFile* out = new TFile("B12_theo_new.root", "recreate");
    TTree* tt = new TTree("T", "branch info");
    int m_num;
    double m_BR;
    int m_np;
    double m_pe[1];
    int m_pn[1];
    tt->Branch("num", &m_num, "num/I");
    tt->Branch("BR", &m_BR, "BR/D");
    tt->Branch("numPhoton", &m_np, "numPhoton/I");
    tt->Branch("photonE", m_pe, "photonE[numPhoton]/D");
    tt->Branch("photonName", m_pn, "photonName[numPhoton]/I");

    m_num = 0;
    m_BR = 0.971908;
    m_np = 0;
    m_pe[0] = 0;
    m_pn[0] = 0;
    tt->Fill();

    m_num = 1;
    m_BR = 0.012296;
    m_np = 1;
    m_pe[0] = 4.438033;
    m_pn[0] = 4440;
    tt->Fill();

    m_num = 2;
    m_BR = 0.014996;
    m_np = 0;
    m_pe[0] = 0;
    m_pn[0] = 0;
    tt->Fill();

    tt->Write();
    hh0->Write();
    hh1->Write();
    hh2->Write();

    out->Close();

}