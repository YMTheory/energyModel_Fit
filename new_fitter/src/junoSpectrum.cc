#include "junoSpectrum.hh"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

using namespace std;

junoSpectrum::junoSpectrum( int nMaxBins,
                            int nMaxBinsData,
                            int nMaxBr,
                            int nMaxGam,
                            double eMin,
                            double eMax,
                            string name)
{

    std::cout << " eMin         = " << eMin << std::endl;
    std::cout << " eMax         = " << eMax << std::endl;
	std::cout << " nMaxBins     = " << nMaxBins << std::endl;
	std::cout << " nMaxBinsData = " << nMaxBinsData << std::endl;
	std::cout << " nMaxBr       = " << nMaxBr << std::endl;
	std::cout << " nMaxGam      = " << nMaxGam << std::endl;
	
    m_eMin      = eMin;
    m_eMax      = eMax;
	m_nBins     = nMaxBins;
    m_nBinsData = nMaxBinsData;
	m_nBranch   = nMaxBr;
	m_nGam      = nMaxGam;
    m_name      = name;

	m_binCenter = new double[nMaxBins];
	
	m_eTru      = new double*[nMaxBr];
	
	m_eTruAlp   = new double [nMaxBr];

	for (int branchIdx = 0; branchIdx < nMaxBr; branchIdx++)
	{
		m_eTru   [branchIdx] = new double[nMaxBins];
	}

	m_eTruGam   = new double*[nMaxBr];
    m_eTruGamStr= new int*[nMaxBr];
	for (int branchIdx = 0; branchIdx < nMaxBr; branchIdx++)
	{
		m_eTruGam[branchIdx] = new double[nMaxGam];
        m_eTruGamStr[branchIdx] = new int[nMaxGam];
	}

    m_eData = new double[nMaxBinsData];
    m_eDataErr = new double[nMaxBinsData];
}

junoSpectrum::~junoSpectrum()
{
	for (int branchIdx = 0; branchIdx < m_nBranch; branchIdx++)
	{
		delete [] m_eTru   [branchIdx];
	}
	for (int branchIdx = 0; branchIdx < m_nBranch; branchIdx++)
	{
		delete [] m_eTruGam[branchIdx];
	}
	delete [] m_binCenter;
	delete [] m_eTru     ;
	delete [] m_eTruGam  ;
	delete [] m_eTruAlp  ;

}

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
        //m_nGamma[branchIdx] = nGamma;
        for(int gamIdx=0; gamIdx<nGamma; gamIdx++) {
            m_eTruGam[branchIdx][gamIdx] = eGamma[gamIdx];
            m_eTruGamStr[branchIdx][gamIdx] = eGammaStr[gamIdx];
        }

        // beta
	    TH1F* electronHist = (TH1F*)ff->Get(Form("hh%d",branchNumber));
		for (int binIdx=0;binIdx!=m_nBins;binIdx++)
		{
			weight = electronHist->Interpolate(m_binCenter[binIdx]);
			m_eTru[branchIdx][binIdx] = branchRatio * weight;
		}
		delete electronHist;
	}

    delete tt;
    delete ff;
}

void junoSpectrum::InitData()
{
    std::cout << " >>> Loading data " << m_name << "spectrum <<< " << std::endl;
    string fileName = "./data/spectrum/data/" + m_name + "_data.root";
    DataHistTree(fileName);
}

void junoSpectrum::DataHistTree(string fileName)
{
    TFile* ff = new TFile(fileName.c_str());
    if(!ff) cout << "No such data file " << fileName << endl;
    TH1D* sigH = (TH1D*)ff->Get(m_name.c_str());
	for (int i=0;i!=m_nBinsData;i++)
	{
		double content = sigH->GetBinContent(i+1);
		double error   = sigH->GetBinError  (i+1);
		m_eData   [i] = content;
		m_eDataErr[i] = error;
	}
	delete sigH;

    delete ff;
}


void junoSpectrum::LoadData()
{
    // Initialization Part
    m_binWidth = (m_eMax - m_eMin) / double(m_nBins);
    m_fitMinBin = int((m_fitMin - m_eMin)/m_binWidth);
    m_fitMinBin = int((m_fitMax - m_eMin)/m_binWidth);

    // theoretical spectrum
    for (int branchIdx=0; branchIdx<m_nBranch; branchIdx++) {
        // beta continuous spectrum
        for (int i=0; i<m_nBins; i++) {
            m_binCenter[i] = m_eMin + m_binWidth*(i+0.5);
            m_eTru[branchIdx][0] = 0;
        }
        // gamma energy
        for (int gamIdx=0; gamIdx<m_nGam; gamIdx++) {
            m_eTruGam[branchIdx][gamIdx] = 0;
        }
        // alpha energy
        m_eTruAlp[branchIdx] = 0;
    }

    // data spectrum
    for (int i=0; i<m_nBinsData; i++) {
        m_eData[i] = 0;
        m_eDataErr[i] = 0;
    }

    // Data Loading ...
    InitTheo();
    InitData();

}














