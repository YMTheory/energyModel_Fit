void text2TH2D()
{
    string filename = "nFe56_all.txt";

    TH2D* hh = new TH2D("nFe56", "", 5000, 0, 5000, 80, 0, 80);
    string line;
    ifstream in; in.open(filename.c_str());
    if(!in) cout << "No such file" << endl;
    double tmp_elec;
    int line_num = 0;
    while(getline(in, line)) {
        int sub_num = 0;
        istringstream ss(line);
        while(ss>>tmp_elec) {
            hh->SetBinContent(line_num+1, sub_num+1, tmp_elec);
            sub_num++;
        }
    line_num++;
    }

    TFile* out = new TFile("nFe56_all.root", "recreate");
    hh->Write();
    out->Close();
}
