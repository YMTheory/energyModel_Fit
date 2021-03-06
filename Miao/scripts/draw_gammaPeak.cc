void draw_gammaPeak()
{
    TFile* fGammaPdf = new TFile("Gamma_Electron1.root", "RECREATE");

    //TGraph* gGammaPdf = new TGraph();
    TH1D* gGammaPdf = new TH1D("h0","",600,0,12);
    string path = "/Users/yumiao/Documents/Works/github/energyModel_Fit/Miao/data/naked_gamma/";
    string name1 = "primary_";  string name2 = ".txt";
    string name[20] = {"511keV", "718keV", "1021keV","1MeV", "2MeV", "4MeV", "4440keV", "6MeV", "8MeV", "Co60", "Cs137", "nH", "K40", "nC12", "Mn54", "Ge68", "nC12", "Tl208", "O16", "nFe56"};
    for(int i=0; i<20; i++)  {

        gGammaPdf->Reset();

        string filename = path+name1+name[i]+name2;
        cout << filename << endl;
    
        gGammaPdf->SetName(("gamma"+name[i]).c_str());

        ifstream in;
        in.open(filename.c_str());
        if(!in) {  cout << "Can not open " << filename << endl; }
        
        string line;
        double tmp_E, tmp_count; int num = 0;
        while(getline(in,line))  {
            istringstream ss(line);
            ss >> tmp_E >> tmp_count;
            //gGammaPdf->SetPoint(num, tmp_E, tmp_count); num++;
            gGammaPdf->SetBinContent(num+1, tmp_count); num++;
        }


        in.close();

        gGammaPdf->Smooth();
        gGammaPdf->Write();
    }



    //hGammaPdf->Delete();
    fGammaPdf->Close();
}
