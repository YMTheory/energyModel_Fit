#include "junoSpectrum.hh"
#include "junoParameters.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

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
    m_eVis      = new double[nMaxBins];
	
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

    m_nPdfBins = 600;
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
    delete [] m_eVis     ;
	delete [] m_eTru     ;
	delete [] m_eTruGam  ;
	delete [] m_eTruAlp  ;
    delete [] m_eData    ;
    delete [] m_eDataErr ;

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
    double weight;
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


void junoSpectrum::ApplyScintillatorNL()
{
    for (int i=0; i<m_nBins; i++) {
        m_eVis[i] = 0;
    }

    //int newBin;
    //int newBinLow, newBinHig;
	//vector<double> eVisGam;

    //for(int branchIdx=0; branchIdx<m_nBranch; branchIdx++) {
	//	eVisGam.push_back(0);
    //    for (int gamIdx=0; gamIdx<m_nGam; gamIdx++) {
	//		if(m_eTruGam[branchIdx][gamIdx]==0) break;  // No more gamma in such branch
    //        //eVisGam[branchIdx] += EvisGamma(to_string(m_eTruGamStr[branchIdx][gamIdx]))*m_eTruGam[branchIdx][gamIdx] * (1+m_gammaScale);

    //    }   
    //}
}

void junoSpectrum::LoadPrmElecDist()
{
    TFile* file = new TFile(junoParameters::gammaPdf_File.c_str(), "read");
    for (int branchIdx=0; branchIdx<m_nBranch; branchIdx++) {
        
        for (int gamIdx=0; gamIdx<m_nGam; gamIdx++) { 
			if(m_eTruGam[branchIdx][gamIdx]==0) break;  // No more gamma in such branch
            string eTru = to_string(m_eTruGamStr[branchIdx][gamIdx]);
            string pdfName = "gamma"+eTru+"keV";
            TH1D* gGammaPdf = (TH1D*)file->Get(pdfName.c_str());
            if(!gGammaPdf) cout << "No Such Pdf : " << pdfName << endl;
            
            int tmp_PdfMaxEtrue;
            double* tmp_pdfEtrue = new double[m_nPdfBins];
            double* tmp_pdfProb = new double[m_nPdfBins];

            for(int i=0; i<gGammaPdf->GetNbinsX(); i++)  {
                tmp_pdfEtrue[i] = gGammaPdf->GetBinCenter(i);
                tmp_pdfProb[i] = gGammaPdf->GetBinContent(i);
                if (tmp_pdfProb[i] == 0) tmp_PdfMaxEtrue = i;
            }

            mapPdfMaxEtrue.insert(pair<int, int> (m_eTruGamStr[branchIdx][gamIdx], tmp_PdfMaxEtrue));
            mapPdfEtrue.insert(pair<int, double*>(m_eTruGamStr[branchIdx][gamIdx], tmp_pdfEtrue));
            mapPdfProb.insert(pair<int, double*> (m_eTruGamStr[branchIdx][gamIdx], tmp_pdfProb));

            delete gGammaPdf;
        }

    }
    
    delete file;
}

double junoSpectrum::EvisGamma(int Etrue)
{
    int gamPdfMaxEtrue   = mapPdfMaxEtrue[Etrue];
    double* gamPdfEtrue  = mapPdfEtrue[Etrue];
    double* gamPdfProb   = mapPdfProb[Etrue];

    double numerator = 0.; double denominator = 0.;
    for(int iBin=1;  iBin<gamPdfMaxEtrue; iBin++) {
        double E1 = gamPdfEtrue[iBin-1];
        double E2 = gamPdfEtrue[iBin];

        double prob1 = gamPdfProb[iBin-1];
        double prob2 = gamPdfProb[iBin];

        double fNL1 = electronQuench::ScintillatorNL(E1) + electronCerenkov::getCerenkovPE(E1);
        double fNL2 = electronQuench::ScintillatorNL(E2) + electronCerenkov::getCerenkovPE(E2);

        numerator   += ( prob1*E1*fNL1 + prob2*E2*fNL2 ) * (E2-E1) /2.;
        denominator += (prob1*E1 + prob2*E2) * (E2-E1)/ 2.;
    } 

    if(denominator ==0) { cout << " >> Error Happens While CalculateGammaNL <<<" << endl; return 0;}
    return numerator/denominator;

}







