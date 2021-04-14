void text2TH2D()
{
    string filename = "K40-nooptical.txt";

    TH2D* hh1 = new TH2D("K40_elec", "", 5000, 0, 5000, 80, 0, 80);
    TH2D* hh2 = new TH2D("K40_posi", "", 5000, 0, 5000, 80, 0, 80);
    
    string line;
    ifstream in; in.open(filename.c_str());
    if(!in) cout << "No such file" << endl;
    double tmp_elec;
    string tmp_num;
    double num;
    int line_num = 0;
    int counta = 0; int countb = 0;
    int count = 0;
    while(getline(in, line)) {
        int sub_num = 0;
        istringstream ss(line);
        //while(ss>>tmp_elec) {
        while(ss>>tmp_num) {
            count++;
            cout << tmp_num << endl;
            if ( tmp_num.find('b') ) {
                cout << " in b chain" << endl;
                tmp_num.erase(tmp_num.end() - 1);
                num = atof(tmp_num.c_str());
                hh1->SetBinContent(line_num+1, sub_num+1, num);
                countb++;
            }
            if ( tmp_num.find('a') ) {
                cout << " in a chain" << endl;
                tmp_num.erase(tmp_num.end() - 1);
                num = atof(tmp_num.c_str());
                hh2->SetBinContent(line_num+1, sub_num+1, num);
                counta++;
            }
            sub_num++;
        }
    line_num++;
    }

    cout << counta << " " << countb << " " << count << endl;

    TFile* out = new TFile("K40_new.root", "recreate");
    hh1->Write();
    hh2->Write();
    out->Close();
}
