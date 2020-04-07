double Integral_BirkLaw(double kB1, double E);
double getCerenkovPE(double E);
double CalculateGammaNL(double *apar, string name);

void draw_fit_gamma()
{
    // draw observed data ...
    TGraphErrors* mTrueGammaNL = new TGraphErrors();
    vector<string> sources;
    vector<double> Etrue;

    double A = 1481.06;
    ifstream in;
    in.open("./data/naked_gamma/singleGamma.txt");
    if(!in){
        cout << " >>> Fail to Open Gamma totPE File !! <<< " << endl;
    }
    string line;

    int count = 0;
    string tmp_name; double tmp_E, tmp_totpe,tmp_totpe_err, tmp_sigma, tmp_sigma_err;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_name >> tmp_E >> tmp_totpe >> tmp_totpe_err  >> tmp_sigma >> tmp_sigma_err;
        sources.push_back(tmp_name);
        Etrue.push_back(tmp_E);
        mTrueGammaNL->SetPoint(count, tmp_E, tmp_totpe/A/tmp_E);
        mTrueGammaNL->SetPointError(count, 0, tmp_sigma/A/tmp_E);
        count++;
    }
    in.close();

    TCanvas* cc = new TCanvas();
    cc->SetGrid();

    mTrueGammaNL->SetMarkerStyle(20);
    mTrueGammaNL->SetLineColor(kBlue+1);
    mTrueGammaNL->SetMarkerColor(kBlue+1);
    mTrueGammaNL->SetLineWidth(2);
    mTrueGammaNL->Draw("APLZ");

    // Calculate Predictions ...
    TGraph* mFitGammaNL = new TGraph();
    //double par[3] = {9.59704e-01, 2.64007e-03, 1.17192e+00};
    //double par[3] = {9.79682e-01, 4.44705e-03,8.79758e-01};
    double par[3] = {0.9798, 6.5e-3, 1.0};
    const int nPoints = sources.size();
    for(int iPoint=0; iPoint<nPoints; iPoint++){
        mFitGammaNL->SetPoint(iPoint, Etrue[iPoint], CalculateGammaNL(par, sources[iPoint])) ;
    }

    mFitGammaNL->SetMarkerStyle(21);
    mFitGammaNL->SetLineColor(kGreen+1);
    mFitGammaNL->SetMarkerColor(kGreen+1);
    mFitGammaNL->SetLineWidth(2);
    mFitGammaNL->Draw("PL SAME");

    TLegend *led = new TLegend();
    led->AddEntry(mTrueGammaNL, "data");
    led->AddEntry(mFitGammaNL, "fitting");
    led->Draw("SAME");

    TLatex latex;
    latex.DrawLatex(3,0.95, "kA=9.79682e-01");
    latex.DrawLatex(3,0.91, "kB=4.44705e-03");
    latex.DrawLatex(3,0.87, "kC=8.79758e-01");
    //latex.Draw("SAME");

    TText text;
    text.DrawText(2.5,1, "Co60");
    text.DrawText(0.662,1, "Cs137");
    text.DrawText(0.8,1, "K40");
    text.DrawText(1.022, 1, "Ge60");
    text.DrawText(1.4, 1, "Mn54");
    text.DrawText(2.2, 1.01, "nH");
    text.Draw("SAME");
}



double Integral_BirkLaw(double kB1, double E)
{
    vector<double> Etrue; vector<double> StopPow;
    ifstream in;
    in.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/gamma/StopPow1.txt");
    if(!in){
        cout << " >>> Fail to Open Quench File!! <<< " << endl;
    }
    string line;

    double tmp_E, tmp_StopPow;
    while(getline(in,line))
    {
        istringstream ss(line);
        ss >> tmp_E >> tmp_StopPow;
        Etrue.push_back(tmp_E);
        StopPow.push_back(tmp_StopPow); 
    }
    in.close();
    
    if(Etrue.size() == 0 or StopPow.size() == 0){
        cout << " >>>  No Data in Vectors! <<< " << endl; return -1;
    } else if(Etrue.size() != StopPow.size() ) {
        cout << " >>> Vector Length are Different! <<< " << endl;  return -1;
    } else {

        // numerical integral to get quenched nonlinearity...
        int num = Etrue.size();
        double sum = 0.;
        double integ_part1 = 0.; double integ_part2 = 0.;
        for(int i=1; i<num; i++){
            integ_part1 = 1./(1+kB1*StopPow[i-1]);
            integ_part2 = 1./(1+kB1*StopPow[i]);
            if( Etrue[i] <= E ){ sum+=(integ_part1+integ_part2)*(Etrue[i]-Etrue[i-1])/2.; }
            else {break;}
        }
        
        return sum/E;
    }
}


double getCerenkovPE(double E)
{
    vector<double> Etrue; vector<double> Cerenkov;
    ifstream in;
    in.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/absolutePE.txt");
    if(!in){
        cout << " >>> Fail to Open Cerenkov File !! <<< " << endl;
    }
    string line;

    double tmp_Edep, tmp_rel, tmp_abs;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_Edep >> tmp_rel >> tmp_abs ;
        Etrue.push_back(tmp_Edep/1000.);
        Cerenkov.push_back(tmp_abs);
    }

    in.close();
    
    if(Cerenkov.size() == 0) {
        cout << " >>> No Data in Cerenkov Vector <<< " << endl; return -1;
    } else if (Cerenkov.size() != Etrue.size()){
        cout << " >>> Vector Length are Different !! <<< " << endl;  return -1;
    } else {

        // get Cerenkov PE
        int num = Cerenkov.size();
        for(int i=1; i<num; i++){
            if(Etrue[i-1]<=E and Etrue[i]>E){  return Cerenkov[i]; }
        }
        cout << E << "   >>> Energy Beyond Range !! <<< " << endl; return -1;
    }
    
}


double CalculateGammaNL(double* apar, string name)
{
    double A = 1481.06;
    double kA = apar[0];
    double kB = apar[1];
    double kC = apar[2];

    // read primary e- dist ...
    
    vector<double> Etrue;
    vector<double> prm_count;

    string name1 = "./data/naked_gamma/primary_";
    string name2 = ".txt";
    string filename = name1 + name + name2;

    ifstream in;
    in.open(filename.c_str());
    if(!in){ cout << "Error Open Primary e+- File" << endl; }
    string line;
    
    double tmp_E, tmp_count;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_E >> tmp_count;
        Etrue.push_back(tmp_E);
        prm_count.push_back(tmp_count);
    }
    in.close();

    const int nBins = Etrue.size();
    double numerator = 0.; double denominator = 0.;
    for(int iBin=1; iBin<nBins-1; iBin++){

        double EE1 = Etrue[iBin-1];
        double EE2 = Etrue[iBin];

        double fq1 = Integral_BirkLaw(kB, EE1) ;
        double fc1 = getCerenkovPE(EE1);
        double electronNL1 = kA*fq1+kC*fc1/A/EE1;
        double fq2 = Integral_BirkLaw(kB, EE2);
        double fc2 = getCerenkovPE(EE2);
        double electronNL2 = kA*fq2+kC*fc2/A/EE2;

        numerator +=   (prm_count[iBin-1]* EE1 *electronNL1 +prm_count[iBin] *EE2*electronNL2 ) * (EE2-EE1) /2.;
        denominator += (prm_count[iBin-1]*EE1 + prm_count[iBin]*EE2)*(EE2-EE1)/2.;
    }

    return numerator/denominator;
}
