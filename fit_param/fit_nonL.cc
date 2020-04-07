void fit_nonL()
{
    ifstream in;
    in.open("../resolution/calib_sources/singleGamma.txt");
    string line;
    const int num = 4;
    double Etrue[num], totPE[num], totPE_err[num], sigma[num], sigma_err[num];

    int count = 0;
    while(getline(in, line)){
       istringstream ss(line);
        ss >> Etrue[count] >> totPE[count] >> totPE_err[count] >> sigma[count] >> sigma_err[count];
        count++;
    }

    TGraphErrors* g1 = new TGraphErrors();

    double energy_scale = 1481.06;   // choose as electron 1MeV energy scale...
    double ratio[num], ratio_err[num];
    for(int i=0; i<num; i++){ 
        ratio[i] = totPE[i]/energy_scale/Etrue[i]; 
        ratio_err[i] = sigma[i]/energy_scale/Etrue[i];
        g1->SetPoint(i, Etrue[i], ratio[i]);
        g1->SetPointError(i, 0, ratio_err[i]);
    }

    TCanvas* cc = new TCanvas();
    cc->SetGrid();

    g1->GetYaxis()->SetRangeUser(0.85,1.05);

    g1->SetMarkerStyle(21);
    g1->SetMarkerColor(kBlue);
    g1->SetLineColor(kBlue);
    g1->SetLineWidth(2);
    //g1->SetTitle("Energy resolution ;E_{vis}/MeV;resolution");
    g1->Draw("AP");


}
