void draw_Quenching()
{
    vector<double> Etrue; vector<double> NonL;
    Etrue.clear();
    NonL.clear();

    TFile* fQuench = new TFile("Quench.root", "RECREATE");
    
    TGraph* gQuench  = new TGraph() ;

    // read stopping power from estar ...
    ifstream in;
    string path = "/Users/yumiao/Documents/Works/Simulation/Nonlinearity/gamma/QuenchFile/kB";
    string back = ".txt";
    for(int i=10; i<100; i++){
        gQuench->Set(0);

        string filename = path+to_string(i)+back;
        cout  << "Open File : " << filename << endl;
        in.open(filename.c_str());
        if(!in){
            cout << " >>> Fail to Open Quench File!! <<< " << endl;
        }
        string line;

        double tmp_E, tmp_NonL; int num = 0;
        while(getline(in,line))
        {
            istringstream ss(line);
            ss >> tmp_E >> tmp_NonL;
            //Etrue.push_back(tmp_E);
            //NonL.push_back(tmp_NonL);
            gQuench->SetPoint(num, tmp_E, tmp_NonL);  num++;
        }
        in.close();
    
        string graph_name = "kB"+to_string(i);
        gQuench->SetName(graph_name.c_str());
        gQuench->SetMarkerStyle(20);
        gQuench->SetMarkerSize(0.5);
        //gQuench->Draw("APL");
        gQuench->Write();

    }

    fQuench->Close();
}
