void test(){
    TGraph* g1 = new TGraph();
    ifstream in;
    in.open("./log");
    string line; double tmp_E, tmp_nl; int i=0;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_E >> tmp_nl;
        g1->SetPoint(i, tmp_E, tmp_nl); i++;
    }
    g1->Draw("AP");

    double p0 = 1.02561e+00;
    double p1 = 1.12245e-01;
    double p2 = 1.39421e+00;
    double p3 = 5.55117e-04;

    TGraph* g2 = new TGraph(); i = 0;
    for(int i=0; i<10000; i++) {
        double xx = 16/10000.*i ;
        double yy =(p0+p3*xx)/(1+p1*TMath::Exp(-p2*xx));
        g2->SetPoint(i, xx, yy);
    }
    
    g2->SetMarkerColor(kRed);
    //g2->Draw("P SAME");
    g2->Draw("AP");
}

