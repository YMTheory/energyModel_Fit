void draw_contour()
{
    gStyle->SetOptStat(0);
    double chimin = 3.16001;
    double m_bestFit[3] = { 9.86065e-01, 6.00267e-03, 8.54511e-01};
 
    // Draw kB, kC ... 
    
    TH2F* contour = new TH2F("contour", "",154, 0.005232, 0.006772, 340, 0.6853, 1.025);
    TGraph* sigma1 = new TGraph(); int num = 0;

    // read chi2 :
    ifstream in;
    in.open("../chi2_p1p2");
    string line;
    double p1, p2, chi2;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> p1 >> p2 >> chi2;
        int binX = int( (p1-0.005232+0.0000001)*100000 );
        int binY = int((p2-0.6853+0.00001)*1000);
        //cout << p1 << " " << binX << " " <<p2 << " " << binY << " " << chi2 << endl;
        contour->SetBinContent(binX, binY, chi2);
        if(TMath::Abs(chi2-chimin)<=0.01)  {
            sigma1->SetPoint(num, p1, p2); num++;
        }
    }

    TGraph* gBest = new TGraph();
    gBest->SetPoint(0, m_bestFit[1], m_bestFit[2]);
    gBest->SetMarkerColor(kRed+1);
    gBest->SetMarkerStyle(29);
    gBest->SetMarkerSize(2);

    sigma1->SetMarkerColor(kRed+1);
    sigma1->SetLineColor(kRed+1);
    sigma1->SetLineWidth(2);
    sigma1->SetMarkerStyle(20);
    sigma1->SetMarkerSize(0.5);

    TText* text1 = new TText(m_bestFit[1]+0.0001, m_bestFit[2], "best fit");
    text1->SetTextColor(kRed+1);
    text1->SetTextFont(43);
    text1->SetTextSize(20);

    TText* text2 = new TText(m_bestFit[1]+0.0001, m_bestFit[2]+0.10, "1 sigma");
    text2->SetTextColor(kRed+1);
    text2->SetTextFont(43);
    text2->SetTextSize(20);

    contour->SetTitle(" ;kB; kC");
    contour->Draw("COLZ");
    gBest->Draw("P SAME");
    sigma1->Draw("P SAME");
    text1->Draw("SAME");
    text2->Draw("SAME");
    

    // Draw kA, kC ...
/*
    TH2F* contour = new TH2F("contour", "", 108, 0.9806, 0.9914, 340, 0.6853, 1.0253);
    TGraph* sigma1 = new TGraph(); int num = 0;

    // read chi2 :
    ifstream in;
    in.open("../chi2_p0p2");
    string line;
    double p1, p2, chi2;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> p1 >> p2 >> chi2;
        int binX = int( (p1-0.9806+0.0000001)*10000 );
        int binY = int((p2-0.6852+0.00001)*1000);
        //cout << p1 << " " << binX << " " <<p2 << " " << binY << " " << chi2 << endl;
        contour->SetBinContent(binX, binY, chi2);

        if(TMath::Abs(chi2-chimin)<=0.01)  {
            sigma1->SetPoint(num, p1, p2); num++;
        }
    }

    TGraph* gBest = new TGraph();
    gBest->SetPoint(0, 9.86065e-01, 8.54511e-01);
    gBest->SetMarkerColor(kRed+1);
    gBest->SetMarkerStyle(29);
    gBest->SetMarkerSize(2);

    sigma1->SetMarkerColor(kRed+1);
    sigma1->SetLineColor(kRed+1);
    sigma1->SetLineWidth(2);
    sigma1->SetMarkerStyle(20);
    sigma1->SetMarkerSize(0.5);

    TText* text1 = new TText(0.9865, 0.8545, "best fit");
    text1->SetTextColor(kRed+1);
    text1->SetTextFont(43);
    text1->SetTextSize(20);

    TText* text2 = new TText(0.9865, 0.955, "1 sigma");
    text2->SetTextColor(kRed+1);
    text2->SetTextFont(43);
    text2->SetTextSize(20);

    contour->SetTitle(" ;kA; kC");
    contour->Draw("COLZ");
    gBest->Draw("P SAME");
    sigma1->Draw("P SAME");
    text1->Draw("SAME");
    text2->Draw("SAME");
*/

    // Draw kA, kB ...
/*
    TH2F* contour = new TH2F("contour", "", 180, 0.9770, 0.9950, 328, 4.01e-3, 7.29e-3);
    TGraph* sigma1 = new TGraph(); int num = 0;

    // read chi2 :
    ifstream in;
    in.open("../chi2_p0p1");
    string line;
    double p1, p2, chi2;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> p1 >> p2 >> chi2;
        int binX = int( (p1-0.9770+0.0000001)*10000 );
        int binY = int((p2-4.01e-3+0.0000001)*100000);
        contour->SetBinContent(binX, binY, chi2);

        if(TMath::Abs(chi2-3.16001)<=0.01)  {
            sigma1->SetPoint(num, p1, p2); num++;
        }
    }

    TGraph* gBest = new TGraph();
    gBest->SetPoint(0, 9.86065e-01, 6.00267e-03);
    gBest->SetMarkerColor(kRed+1);
    gBest->SetMarkerStyle(29);
    gBest->SetMarkerSize(2);

    sigma1->SetMarkerColor(kRed+1);
    sigma1->SetLineColor(kRed+1);
    sigma1->SetLineWidth(2);
    sigma1->SetMarkerStyle(20);
    sigma1->SetMarkerSize(0.5);

    TText* text1 = new TText(0.9865, 6.00e-3, "best fit");
    text1->SetTextColor(kRed+1);
    text1->SetTextFont(43);
    text1->SetTextSize(20);

    TText* text2 = new TText(0.9865, 6.64e-3, "1 sigma");
    text2->SetTextColor(kRed+1);
    text2->SetTextFont(43);
    text2->SetTextSize(20);

    contour->SetTitle(" ;kA; kB");
    contour->Draw("COLZ");
    gBest->Draw("P SAME");
    sigma1->Draw("P SAME");
    text1->Draw("SAME");
    text2->Draw("SAME");
*/
}
