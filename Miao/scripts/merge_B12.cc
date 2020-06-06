void merge_B12()
{
    TH1D* hh0 = new TH1D("hh0", "", 100, 0, 14);
    TH1D* hh1 = new TH1D("hh1", "", 100, 0, 14);
    TH1D* hh2 = new TH1D("hh2", "", 100, 0, 14);

    double bin0[100] = {62.3567,193.02,344.337,515.528,705.151,911.679,1133.59,1369.37,1617.56,1876.71,2145.43,2422.35,2706.15,2995.52,3289.23,3586.06,3884.82,4184.39,4483.66,4781.57,5077.09,5369.25,5657.08,5939.69,6216.19,6485.77,6747.61,7000.97,7245.13,7479.4,7703.15,7915.77,8116.69,8305.4,8481.4,8644.24,8793.51,8928.84,9049.88,9156.36,9248,9324.59,9385.95,9431.93,9462.42,9477.37,9476.74,9460.55,9428.83,9381.68,9319.22,9241.61,9149.05,9041.79,8920.1,8784.28,8634.7,8471.75,8295.84,8107.46,7907.09,7695.29,7472.62,7239.71,6997.21,6745.82,6486.24,6219.27,5945.68,5666.34,5382.1,5093.9,4802.67,4509.4,4215.12,3920.89,3627.8,3337,3049.64,2766.94,2490.13,2220.5,1959.37,1708.07,1468,1240.58,1027.26,829.543,648.949,487.045,345.426,225.728,129.616,58.7937,14.9979,0,0,0,0,0};
    double bin1[100] = {27.7951,85.1174,150.176,222.301,300.544,383.939,471.542,562.441,655.758,750.655,846.329,942.014,1036.98,1130.55,1222.07,1310.92,1396.53,1478.36,1555.93,1628.77,1696.45,1758.61,1814.89,1864.99,1908.64,1945.62,1975.73,1998.83,2014.79,2023.55,2025.06,2019.32,2006.39,1986.32,1959.24,1925.3,1884.7,1837.65,1784.43,1725.35,1660.74,1590.99,1516.52,1437.78,1355.26,1269.5,1181.07,1090.57,998.644,905.979,813.29,721.334,630.903,542.828,457.975,377.249,301.589,231.975,169.419,114.973,69.7239,34.7956,11.348,0.577355,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double bin2[100] = {11.2863,33.9298,58.7228,85.1978,112.795,140.969,169.209,197.036,224.01,249.73,273.83,295.986,315.908,333.349,348.096,359.977,368.857,374.639,377.265,376.716,373.009,366.2,356.385,343.696,328.302,310.414,290.277,268.175,244.431,219.406,193.497,167.139,140.807,115.011,90.3,67.2596,46.5135,28.7227,14.5853,4.83688,0.250074,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    double scale1 = 0; double scale2 = 0; double scale0 = 0;
    for(int i=0; i<100; i++) {
       scale0 += bin0[i]; scale1 += bin1[i]; scale2 += bin2[i]; 
    }

    for(int i=0; i<100; i++) {
        hh0->SetBinContent(i+1, bin0[i]);
        hh1->SetBinContent(i+1, bin1[i]/scale1*scale0);
        hh2->SetBinContent(i+1, bin2[i]/scale2*scale0);
    }

    double b0 = 0.98;
    double b1 = 0.01182;
    double b2 = 1-b0-b2;


    TTree * newTree = new TTree("T","T");
    Int_t m_num;
    newTree->Branch("num", &m_num, "num/I");
    Double_t m_branch;
    newTree->Branch("BR", &m_branch, "BR/D");
    Int_t m_numPhoton;
    newTree->Branch("numPhoton", &m_numPhoton, "numPhoton/I");
    Double_t m_photonE[100];
    newTree->Branch("photonE", &m_photonE, "photonE[numPhoton]/D");
    int m_photonName[100];
    newTree->Branch("photonName", m_photonName, "photonName[numPhoton]/I");

    m_num = 0;
    m_branch = b0;
    m_numPhoton = 0;
    newTree->Fill();

    m_num = 1;
    m_branch = 0.01182;
    m_numPhoton = 1;
    m_photonE[0] = 4.44;
    m_photonName[0] = 4440;
    newTree->Fill();

    m_num = 2;
    m_branch = 0.00518;
    m_numPhoton = 0;
    newTree->Fill();


    TFile* file = new TFile("B12_theo.root", "recreate");
    hh0->Write();
    hh1->Write();
    hh2->Write();
    newTree->Write();
    file->Close();


    //TH1D* hh = new TH1D("B12", "", 100, 0, 13.4);
    //for(int i=0; i<100; i++) {
    //    double c0 = hh0->GetBinContent(i+1);
    //    double c1 = hh1->GetBinContent(i+1);
    //    double c2 = hh2->GetBinContent(i+1);
    //    hh->AddBinContent(i+1, b0*c0);
    //    hh->AddBinContent(i+1, b1*c1);
    //    hh->AddBinContent(i+1, b2*c2);
    //    scale1 += hh->GetBinContent(i+1);
    //}

    //hh->SetStats(0);
    //hh->SetLineColor(kBlue);
    //hh->SetLineWidth(2);
    //hh->SetMarkerSize(0.8);
    //hh->SetMarkerColor(kBlue+1);
    //hh->SetMarkerStyle(20);
    //hh->SetTitle("B12 Spectrum; beta kinetic energy/MeV; a.u");
    //hh->Draw("PEX0");
    
    //TH1D* hSimul = new TH1D("Simul", "", 100, 0, 13.4);
    //ifstream in;
    //in.open("./B12_edep.txt");
    //string line; double tmp;
    //while(getline(in, line)) {
    //    istringstream ss(line);
    //    ss >> tmp;
    //    hSimul->Fill(tmp);
    //}
    //double scale2 = hSimul->GetEntries();
    //hSimul->Scale(scale1/scale2);
    //hSimul->SetLineWidth(2);
    //hSimul->SetLineColor(kMagenta);
    //hSimul->SetMarkerColor(kMagenta);
    //hSimul->SetMarkerStyle(20);
    //hSimul->SetMarkerSize(0.8);
    //hSimul->Draw("PEX0 SAME");

    //TLegend* led = new TLegend();
    //led->AddEntry(hh2, "Calculation", "L");
    //led->AddEntry(hSimul, "gendecay sim", "L");
    //led->Draw("SAME");
}
