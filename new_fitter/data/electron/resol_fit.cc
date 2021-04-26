void resol_fit()
{
    TGraphErrors* g1 = new TGraphErrors();

    ifstream in; in.open("./elecResol1.txt");
    string line;
    double Etrue, totpe, totpe_err, sigma, sigma_err, resol, resol_err;
    int index = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> Etrue >> totpe >> totpe_err >> sigma >> sigma_err >> resol >> resol_err ;
        g1->SetPoint(index, Etrue, sigma*sigma);
        g1->SetPointError(index, 0, sigma_err*2*sigma_err);
        index++;
    }

    TF1* f1 = new TF1("f1", "[0] + [1]*x + [2]*x*x", 0, 8);
    g1->Fit(f1, "RE");
    g1->SetMarkerSize(0.45);
    g1->SetMarkerColor(kBlue+1);
    g1->SetLineColor(kBlue+1);
    
    g1->Draw("AP");
}
