void test_gaus()
{
    gStyle->SetOptFit(1111);

    TH1D* h0 = new TH1D("h0", "", 100, 0, 2000);
    TH1D* h1 = new TH1D("h1", "", 100, 0, 2000);
    TH1D* h2 = new TH1D("h2", "", 100, 0, 100);

    Double_t mean[101];
    Double_t sigma[101];
    for(int i=0; i<101; i++) {
        mean[i] = gRandom->Gaus(1000,50);
        sigma[i] = gRandom->Gaus(50,4);
        h0->Fill(mean[i]);
        h2->Fill(sigma[i]);
    }
    
    for(int i=0; i<10000; i++) {
        int idx = int(gRandom->Uniform(0, 101));
        double tmp = gRandom->Gaus(mean[idx], sigma[idx]);
        h1->Fill(tmp);
    }


    TCanvas* c1 = new TCanvas("c1");
    h1->Fit("gaus");
    h1->Draw();

    TCanvas* c2 = new TCanvas("c2");
    h0->Draw("");

    TCanvas* c3 = new TCanvas("c3");
    h2->Draw("");
}
