void toyFit()
{
    TH1D* h1 = new TH1D("h1", "", 100, 900, 1100);
    TH1D* h2 = new TH1D("h2", "", 200, 900, 1100);

    for (int i=0; i<10000; i++) {
        double sp = gRandom->Gaus(1000, 10);
        h1->Fill(sp);
        h2->Fill(sp);
    }

    TF1* f1 = new TF1("f1", "[0]/TMath::Sqrt(2*TMath::Pi())/[2] * TMath::Exp(-(x-[1])*(x-[1])/2/[2]/[2])", 900, 1100);
    TF1* f2 = new TF1("f2", "[0]/TMath::Sqrt(2*TMath::Pi())/[2] * TMath::Exp(-(x-[1])*(x-[1])/2/[2]/[2])", 900, 1100);

    h1->Fit("gaus");
    h2->Fit("gaus");
    //    
    //cout << f1->GetParameter(0) << " " << f1->GetParameter(1) << " " << f1->GetParameter(2) << endl;
    //cout << f2->GetParameter(0) << " " << f2->GetParameter(1) << " " << f2->GetParameter(2) << endl;


    h1->SetLineColor(kBlue+1);
    h1->Draw();
    h2->SetLineColor(kRed+1);
    h2->Draw("same");


}
