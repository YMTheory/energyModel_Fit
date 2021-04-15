void prmelec_check()
{

    string filelist[8];
    filelist[0] = "./Cs137_new.root";
    filelist[1] = "./Mn54_new.root";
    filelist[2] = "./K40_new.root";
    filelist[3] = "./nH_new.root";
    filelist[4] = "./Tl208_new.root";
    filelist[5] = "./nC12_new.root";
    filelist[6] = "./O16_new.root";
    filelist[7] = "./nFe56_new.root";

    string source[8] = {"Cs137", "Mn54", "K40", "nH", "Tl208", "nC12", "O16", "nFe56"};

    double Etrue[8] = {0.662, 0.834, 1.461, 2.223, 2.614, 4.945, 6.130, 7.637};

    TGraphErrors* g1 = new TGraphErrors();

    //TH1D *h1 = new TH1D("edep", "", 10000, 0, 8);
    TH1D *h1 = new TH1D("edep", "", 50, 0, 50);

    for (int i = 0; i < 8; i++)
    {
        h1->Reset();

        string filename = filelist[i];

        TFile *ff = new TFile(filename.c_str(), "read");
        TH2D *hh1 = (TH2D *)ff->Get((source[i] + "_elec").c_str());
        TH2D *hh2 = (TH2D *)ff->Get((source[i] + "_posi").c_str());

        for (int binx = 0; binx < hh1->GetNbinsX(); binx++)
        {
            double evtE = 0;
            double evtN = 0;
            for (int biny = 0; biny < hh1->GetNbinsY(); biny++)
            {
                double cont1 = hh1->GetBinContent(binx, biny);
                if (cont1 != 0) evtN++;
                double cont2 = hh2->GetBinContent(binx, biny);
                if (cont2 != 0) { cont2 += 2*0.511; evtN++; }
                evtE += cont1;
                evtE += cont2;
            }
            //h1->Fill(evtE);
            h1->Fill(evtN);
        }

        double edep = h1->GetMean();
        double esigma = h1->GetStdDev();

        //g1->SetPoint(i, Etrue[i], edep/Etrue[i]);
        g1->SetPoint(i, Etrue[i], edep);

        delete hh1;
        delete hh2;
        delete ff;
    }

    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kBlue+1);
    g1->SetLineColor(kBlue+1);
    g1->SetTitle("; Etrue/MeV; total primary e+- number");
    g1->Draw("APL");
}
