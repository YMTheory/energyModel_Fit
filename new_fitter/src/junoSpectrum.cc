#include "junoSpectrum.hh"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

using namespace std;

void junoSpectrum::InitTheo()
{
    std::cout << " >>> Loading theoretical " << m_name << "spectrum <<< " << std::endl;
    string theofile = "./data/spectrum/theo/" + m_name + "_theo.root";
    TheoHistTree(theofile);
}

void junoSpectrum::TheoHistTree(string theofile)
{
    double energyScale = 1.0;

    TFile* ff = new TFile(theofile.c_str());
    if(!ff) cout << " >>> No theoretical spectrum file " << theofile << endl;
    TTree* tt = (TTree*)ff->Get("T");
    double eGamma[10];
    int eGammaStr[10];
    int nGamma, branchNumber;
    double branchRatio;
    double weight, weight1;
	tt->SetBranchAddress("num"      ,&branchNumber);
	tt->SetBranchAddress("BR"       ,&branchRatio);
	tt->SetBranchAddress("numPhoton",&nGamma);
	tt->SetBranchAddress("photonE"  ,eGamma);
    tt->SetBranchAddress("photonName", eGammaStr);
    m_nBranch = tt->GetEntries();
    cout << " >>> Total Branch Number = " << m_nBranch << endl;
	for (int branchIdx=0;branchIdx!=tt->GetEntries();branchIdx++){ 
        tt->GetEntry(branchIdx);
        cout << " >>> " << branchIdx << " with "<< nGamma << " gamma and branch ratio is " << branchRatio << endl;
        // gamma from this branch :
        m_nGamma[branchIdx] = nGamma;
        for(int gamIdx=0; gamIdx<nGamma; gamIdx++) {
            m_eTruGam[branchIdx][gamIdx] = eGamma[gamIdx];
            m_eTruGamStr[branchIdx][gamIdx] = eGammaStr[gamIdx];
        }

        // beta
	    TH1F* electronHist = (TH1F*)ff->Get(Form("hh%d",branchNumber));
		for (int binIdx=0;binIdx!=m_nBins;binIdx++)
		{
			weight = electronHist->Interpolate(m_binCenter[binIdx]*energyScale);
			if(isotope==0) 
				m_eTru[branchIdx][binIdx] = branchRatio * weight;
			if(isotope==1)
				m_eTruBck[branchIdx][binIdx] = branchRatio * weight;
		}
		delete electronHist;
	}

    delete tt;
    delete ff;
}
