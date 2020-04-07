double Integral_BirkLaw(double kB1, double E);
double getCerenkovPE(double E);

void draw_fit_electron()
{
    // draw observed data ...
    TGraphErrors* mTrueElectronNL = new TGraphErrors();

    double A = 1481.06;
    ifstream in;
    in.open("./data/electron_total.txt");
    if(!in){
        cout << " >>> Fail to Open Electron totPE File !! <<< " << endl;
    }
    string line;

    int count = 0;
    double tmp_E, tmp_totpe, tmp_sigma;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_E >> tmp_totpe  >> tmp_sigma;
        mTrueElectronNL->SetPoint(count, tmp_E, tmp_totpe/A/tmp_E);
        mTrueElectronNL->SetPointError(count, 0, tmp_sigma/A/tmp_E);
        count++;
    }
    in.close();


    // draw fitting ...
    TGraph* mFitElectronNL = new TGraph();
    double kA = 9.799e-01;
    double kB1 = 6.60101e-03;
    double kC = 9.98625e-01;
    //double kA = 0.9796;//9.56497e-01;
    //double kB1 = 6.5e-03; //6.51875e-06;
    //double kC = 1.e+00;
    double pred_nl[798];
    for(int i=0; i<798; i++){
        double energy = 7.98/798*(i+1);
        //pred_nl[i] =  getCerenkovPE(energy);
        pred_nl[i] = kA*Integral_BirkLaw(kB1,  energy) + kC*getCerenkovPE(energy)/A/energy;
        mFitElectronNL->SetPoint(i, energy, pred_nl[i]);
    }

    //mFitElectronNL->GetYaxis()->SetRangeUser(0.95,1.05);

    TCanvas* cc = new TCanvas();
    cc->SetGrid();
    //gPad->SetLogx();

    mFitElectronNL->SetTitle("Electron Nonlinearity Fitting; Etrue/MeV; Evis/Etrue");

    mTrueElectronNL->SetMarkerStyle(24);
    mTrueElectronNL->SetMarkerSize(0.8);
    mTrueElectronNL->SetMarkerColor(kBlue+1);
    mTrueElectronNL->SetLineColor(kBlue+1);
    mTrueElectronNL->SetLineWidth(2);

    mFitElectronNL->SetMarkerStyle(20);
    mFitElectronNL->SetMarkerSize(0.3);
    mFitElectronNL->SetMarkerColor(kGreen+1);
    mFitElectronNL->SetLineColor(kGreen+1);
    mFitElectronNL->SetLineWidth(3);
    mFitElectronNL->Draw("AL");
    mTrueElectronNL->Draw("PZ SAME");

    TLegend* le = new TLegend(.75,.80,.95,.95);
    le->AddEntry(mTrueElectronNL, "data");
    le->AddEntry(mFitElectronNL, "fitting");
    le->Draw("SAME");

    TLatex latex;
    latex.DrawLatex(3,0.95, "kA=9.79991e-01");
    latex.DrawLatex(3,0.91, "kB=6.60101e-03");
    latex.DrawLatex(3,0.87, "kC=9.98625e-01");
    latex.Draw("SAME");

    return;
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
