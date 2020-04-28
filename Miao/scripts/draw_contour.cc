void draw_contour()
{
    gStyle->SetOptStat(0);
 
    // Draw kB, kC ... 
   /* 
    TH2F* contour = new TH2F("contour", "", 400, 4e-3, 8e-3, 800, 1.01, 1.09);
    TGraph* sigma1 = new TGraph(); int num = 0;

    // read chi2 :
    ifstream in;
    in.open("../chi2_p1p2");
    string line;
    double p1, p2, chi2;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> p1 >> p2 >> chi2;
        int binX = int( (p1-4e-3+0.0000001)*100000 );
        int binY = int((p2-1.01+0.00001)*10000);
        //cout << p1 << " " << binX << " " <<p2 << " " << binY << " " << chi2 << endl;
        contour->SetBinContent(binX, binY, chi2);
        if(TMath::Abs(chi2-1.356483)<=0.01)  {
            sigma1->SetPoint(num, p1, p2); num++;
        }
    }

    TGraph* gBest = new TGraph();
    gBest->SetPoint(0, 5.11806e-03, 1.04938e+00);
    gBest->SetMarkerColor(kRed+1);
    gBest->SetMarkerStyle(29);
    gBest->SetMarkerSize(2);

    sigma1->SetMarkerColor(kRed+1);
    sigma1->SetLineColor(kRed+1);
    sigma1->SetLineWidth(2);
    sigma1->SetMarkerStyle(20);
    sigma1->SetMarkerSize(0.5);

    TText* text1 = new TText(5.118e-3, 1.043, "best fit");
    text1->SetTextColor(kRed+1);
    text1->SetTextFont(43);
    text1->SetTextSize(20);

    TText* text2 = new TText(5.2e-3, 1.073, "1 sigma");
    text2->SetTextColor(kRed+1);
    text2->SetTextFont(43);
    text2->SetTextSize(20);

    contour->SetTitle(" ;kB; kC");
    contour->Draw("COLZ");
    gBest->Draw("P SAME");
    sigma1->Draw("P SAME");
    text1->Draw("SAME");
    text2->Draw("SAME");
 */   

    // Draw kA, kC ...
/*
    TH2F* contour = new TH2F("contour", "", 600, 0.970, 0.976, 800, 1.01, 1.09);
    TGraph* sigma1 = new TGraph(); int num = 0;

    // read chi2 :
    ifstream in;
    in.open("../chi2_p0p2");
    string line;
    double p1, p2, chi2;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> p1 >> p2 >> chi2;
        int binX = int( (p1-0.970+0.0000001)*100000 );
        int binY = int((p2-1.01+0.00001)*10000);
        //cout << p1 << " " << binX << " " <<p2 << " " << binY << " " << chi2 << endl;
        contour->SetBinContent(binX, binY, chi2);

        if(TMath::Abs(chi2-1.356483)<=0.001)  {
            sigma1->SetPoint(num, p1, p2); num++;
        }
    }

    TGraph* gBest = new TGraph();
    gBest->SetPoint(0, 0.973, 1.04938e+00);
    gBest->SetMarkerColor(kRed+1);
    gBest->SetMarkerStyle(29);
    gBest->SetMarkerSize(2);

    sigma1->SetMarkerColor(kRed+1);
    sigma1->SetLineColor(kRed+1);
    sigma1->SetLineWidth(2);
    sigma1->SetMarkerStyle(20);
    sigma1->SetMarkerSize(0.5);

    TText* text1 = new TText(0.9732, 1.049, "best fit");
    text1->SetTextColor(kRed+1);
    text1->SetTextFont(43);
    text1->SetTextSize(20);

    TText* text2 = new TText(0.9732, 1.071, "1 sigma");
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

    TH2F* contour = new TH2F("contour", "", 600, 0.970, 0.976, 400, 4e-3, 8e-3);
    TGraph* sigma1 = new TGraph(); int num = 0;

    // read chi2 :
    ifstream in;
    in.open("../chi2_p0p1");
    string line;
    double p1, p2, chi2;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> p1 >> p2 >> chi2;
        int binX = int( (p1-0.970+0.0000001)*100000 );
        int binY = int((p2-4e-3+0.0000001)*100000);
        contour->SetBinContent(binX, binY, chi2);

        if(TMath::Abs(chi2-1.356483)<=0.01)  {
            sigma1->SetPoint(num, p1, p2); num++;
        }
    }

    TGraph* gBest = new TGraph();
    gBest->SetPoint(0, 0.973, 5.11806e-3);
    gBest->SetMarkerColor(kRed+1);
    gBest->SetMarkerStyle(29);
    gBest->SetMarkerSize(2);

    sigma1->SetMarkerColor(kRed+1);
    sigma1->SetLineColor(kRed+1);
    sigma1->SetLineWidth(2);
    sigma1->SetMarkerStyle(20);
    sigma1->SetMarkerSize(0.5);

    TText* text1 = new TText(0.9732, 5.12e-3, "best fit");
    text1->SetTextColor(kRed+1);
    text1->SetTextFont(43);
    text1->SetTextSize(20);

    TText* text2 = new TText(0.9732, 4.64e-3, "1 sigma");
    text2->SetTextColor(kRed+1);
    text2->SetTextFont(43);
    text2->SetTextSize(20);

    contour->SetTitle(" ;kA; kB");
    contour->Draw("COLZ");
    gBest->Draw("P SAME");
    sigma1->Draw("P SAME");
    text1->Draw("SAME");
    text2->Draw("SAME");

}
