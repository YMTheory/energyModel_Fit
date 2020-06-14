void generateFile()
{

    double scale = 1481.06;

    // Read Electron Nonlinearity Curve
    double elecEtrue[799]; double elecNonl[799];
    ifstream in; in.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/MC_NonL.txt");
    string line; double tmp_nonl; int num = 0;
    while(getline(in,line)) {
        istringstream ss(line);
        ss >> tmp_nonl;
        elecEtrue[num] = num*0.01; elecNonl[num] = tmp_nonl; num++;
    }
    in.close();

    // Read Electron Resolution Curve
    double elecTotPE[799]; double elecTotPESigma[799];
    in.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/electron/Cerenkov/totPE_smearing.txt");
    if(!in) {cout << "No such file !!" << endl;} 
    num = 0; double tmp_E, tmp_totpe, tmp_sigma;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_E >> tmp_totpe >> tmp_sigma;
        elecTotPE[num] = tmp_totpe; elecTotPESigma[num] = tmp_sigma; num++;
    }
    in.close();


    vector<double> EprmElec;
    in.open("/Users/yumiao/Documents/Works/Simulation/Nonlinearity/gamma/rootfiles/Cs137.txt");
    if(!in) {cout << "No such file !!" << endl;}
    Double_t tmp_elec; Int_t index = 0;
    Double_t Evis_calc[1000]; 
    Double_t totpe_calc[1000]; Double_t totpeSigma_calc[1000];
    while(getline(in,line)) {
        cout << "Processing " << index << " event" << endl;
        num = 0;
        istringstream ss(line);
        while(ss>>tmp_elec) {
            EprmElec.push_back(tmp_elec); num++;
        }
        double tmp_Evis = 0; double tmp_sigma = 0;
        for(int j=0; j<num; j++) {
            int idx = int(EprmElec[j]/0.01);
            double elec_nonl = elecNonl[idx];
            tmp_Evis += elec_nonl*EprmElec[j];
            tmp_sigma += elecTotPESigma[idx]*elecTotPESigma[idx];
            //tmp_Evis += EprmElec[j];
        }

        Evis_calc[index] = tmp_Evis*scale; 
        totpe_calc[index] = tmp_Evis*scale;
        totpeSigma_calc[index] = TMath::Sqrt(tmp_sigma);
        index++;
        EprmElec.clear();

        if (index==1000) break;
    }
    in.close();

    ofstream out("singleEvent_Cs137.txt", ios::app);
    for(int i=0; i<1000; i++) {
        out << totpe_calc[i] <<" "<<totpeSigma_calc[i] << endl;
    }
    out.close();
}
