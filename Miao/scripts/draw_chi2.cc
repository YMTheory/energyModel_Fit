void draw_chi2()
{
    TGraph* g1 = new TGraph();

    ifstream in; in.open("iteration.txt");
    string line; int i = 0;
    double time, chi2;
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> time >> chi2;
        g1->SetPoint(i, time, chi2); i++;
    }
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kBlue+1);
    g1->SetLineColor(kBlue+1);
    g1->SetMarkerSize(0.5);
    g1->Draw("APL");
}
