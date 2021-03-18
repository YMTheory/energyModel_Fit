#include "junoSpectrum.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "junoParameters.hh"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

#include <iostream>

using namespace std;

double junoSpectrum::m_pdf_eTru[2000];
double junoSpectrum::m_pdf_prob[2000];
double junoSpectrum::m_max_eTru;

double junoSpectrum::m_gammaScale = 0;

junoSpectrum::junoSpectrum( int nMaxBins,
                            int nMaxBinsData,
                            int nMaxBr,
                            int nMaxGam,
                            double eMin,
                            double eMax,
                            double fitMin,
                            double fitMax,
                            string name)

{
	std::cout << " nMaxBins     = " << nMaxBins << std::endl;
	std::cout << " nMaxBinsData = " << nMaxBinsData << std::endl;
	std::cout << " nMaxBr       = " << nMaxBr << std::endl;
	std::cout << " nMaxGam      = " << nMaxGam << std::endl;
	
	m_nBins     = nMaxBins;
	m_nBinsData = nMaxBinsData;
	m_nBr       = nMaxBr;
	m_nGam      = nMaxGam;
    m_eMin      = eMin;
    m_eMax      = eMax;
    m_fitMin    = fitMin;
    m_fitMax    = fitMax;
    m_name = name;
    m_title = name;
	
	m_binCenter = new double[nMaxBins];
	m_eVis      = new double[nMaxBins];
	m_eVisBck   = new double[nMaxBins];
	m_eRec      = new double[nMaxBins];
	m_eRecBck   = new double[nMaxBins];
	m_eTheo     = new double[nMaxBinsData];
	m_eData     = new double[nMaxBinsData];
	m_eDataErr  = new double[nMaxBinsData];
	
	m_eTru      = new double*[nMaxBr];
	m_eTruBck   = new double*[nMaxBr];
	
	m_eTruAlp   = new double [nMaxBr];
	
	for (int branchIdx = 0; branchIdx < nMaxBr; branchIdx++)
	{
		m_eTru   [branchIdx] = new double[nMaxBins];
		m_eTruBck[branchIdx] = new double[nMaxBins];
	}
	m_eTruGam   = new double*[nMaxBr];
    m_eTruGamStr= new int*[nMaxBr];
	for (int branchIdx = 0; branchIdx < nMaxBr; branchIdx++)
	{
		m_eTruGam[branchIdx] = new double[nMaxGam];
        m_eTruGamStr[branchIdx] = new int[nMaxGam];
	}

}


junoSpectrum::~junoSpectrum()
{
	std::cout << " Destroying " << m_name << " <------" << std::endl;
	//std::cout << " **** m_nBr = " << m_nBr << std::endl;
	for (int branchIdx = 0; branchIdx < m_nBr; branchIdx++)
	{
		delete [] m_eTru   [branchIdx];
		delete [] m_eTruBck[branchIdx];
	}
	for (int branchIdx = 0; branchIdx < m_nBr; branchIdx++)
	{
		delete [] m_eTruGam[branchIdx];
	}
	delete [] m_binCenter;
	delete [] m_eVis     ;
	delete [] m_eVisBck  ;
	delete [] m_eRec     ;
	delete [] m_eRecBck  ;
	delete [] m_eTheo    ;
	delete [] m_eData    ;
	delete [] m_eDataErr ;
	delete [] m_eTru     ;
	delete [] m_eTruBck  ;
	delete [] m_eTruGam  ;
	delete [] m_eTruAlp  ;
}

void junoSpectrum::LoadData(string fileName)
{
    m_opt           = true;
    m_dataIsLoaded  = false;

	m_binWidth  = (m_eMax-m_eMin)/double(m_nBins);
	m_fitMinBin = int((m_fitMin-m_eMin)/m_binWidth);
	m_fitMaxBin = int((m_fitMax-m_eMin)/m_binWidth);

	for(int branchIdx=0; branchIdx<m_nBr; branchIdx++)
	{
		for(int i=0; i<m_nBins; i++)
		{
            m_binCenter[i] = m_eMin + m_binWidth*(i+0.5); //cout << "binCenter: " << m_binCenter[i] << endl;
			m_eTru   [branchIdx][i] = 0;
			m_eTruBck[branchIdx][i] = 0;
		}
		for(int gamIdx=0; gamIdx<m_nGam; gamIdx++)
		{
			m_eTruGam[branchIdx][gamIdx] = 0;
		}
		m_eTruAlp[branchIdx] = 0;
	}
    InitTheo();
	//if (fileName.find(".root") != std::string::npos)
	InitData(fileName);
	//else
	//	InitToyMC(fileName);
	m_dataIsLoaded = true;
	cout << " finished preparation of " << m_name << " data <------ " << endl;  
}


void junoSpectrum::TheoHistTree(string filename, int isotope)
{
    double energyScale = 1.0;

	TFile *file = new TFile(filename.c_str());
	TTree *tree = (TTree*)file->Get("T");
    double eGamma[10];
    int eGammaStr[10];
	int nGamma,branchNumber;
	double branchRatio;
	double weight,weight1;
	tree->SetBranchAddress("num"      ,&branchNumber);
	tree->SetBranchAddress("BR"       ,&branchRatio);
	tree->SetBranchAddress("numPhoton",&nGamma);
	tree->SetBranchAddress("photonE"  ,eGamma);
    tree->SetBranchAddress("photonName", eGammaStr);
    //if(m_name=="B12")
    std::cout << " R*******eading theory hist " << std::endl;
    std::cout << "total branch number : " << tree->GetEntries() << endl;
	for (int branchIdx=0;branchIdx!=tree->GetEntries();branchIdx++){
		tree->GetEntry(branchIdx);
		//if(m_name=="B12")
            std::cout << " ------> " << branchIdx  << " with " << nGamma  << " gamma and branch ratio: "<< branchRatio  <<  std::endl;
		/// gammas
		for (int gamIdx=0;gamIdx!=nGamma;gamIdx++)
		{
			m_eTruGam[branchIdx][gamIdx]    = eGamma[gamIdx]/energyScale;
            m_eTruGamStr[branchIdx][gamIdx] = eGammaStr[gamIdx];
			std::cout <<  eGamma[gamIdx] << " " << eGammaStr[gamIdx] << std::endl;
			//if(m_name=="B12")
		}
		/// electrons
		TH1F* electronHist = (TH1F*)file->Get(Form("hh%d",branchNumber));
		
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
    delete tree;
    file->Close();
    delete file;

}

void junoSpectrum::InitData(string fileName)
{
    std::cout << " >>> Loading spectrum data from " << fileName << endl;
    DataHist(fileName);
}

void junoSpectrum::DataHist(string fileName)
{
	TFile* file = new TFile(fileName.c_str());
	TH1F* sigH  = (TH1F*)file->Get("Michel");
    for (int i=0; i!=m_nBinsData;i++)
	{
		double content = sigH->GetBinContent(i+1);
		double error   = sigH->GetBinError  (i+1);
		m_eData   [i] = content;
		m_eDataErr[i] = error;
	}
	delete sigH;
	file->Close();
	delete file;

}


void junoSpectrum::Normalize()
{
	int   rebin = m_nBins/m_nBinsData;
	double binWidthData = m_binWidth * rebin;
	double nTheo = 0;
	double nData = 0;
	for (int i = 0; i < m_nBinsData; i++)
	{
		m_eTheo[i] = 0;
        for (int j = 0; j < rebin; j++){
			m_eTheo[i] += m_eVis[i*rebin+j];
        } 
            
		if(i*binWidthData>m_fitMin && i*binWidthData<m_fitMax)
		{
			nTheo += m_eTheo[i];
			nData += m_eData[i];
		}
	}
    double scale = 1;
    if( nTheo!=0 ) { scale = nData/nTheo; }
	for (int i = 0; i < m_nBinsData; i++)
    {
		m_eTheo[i] *= scale;
        
    }
	for (int i = 0; i < m_nBins; i++)
	{
		m_eVis   [i] *= scale;
		m_eVisBck[i] *= scale;
	}
}

void junoSpectrum::ApplyScintillatorNL()
{
    for(int i=0; i<m_nBins; i++) {
        m_eVis   [i] = 0;
        m_eVisBck[i] = 0;
    }

    int newBin, newBinBck;
	int    newBinLow,newBinLowBck;
	int    newBinHig,newBinHigBck;
	double bias,biasBck;
	double eTru;
	double eVis,eVisBck;
	double eVisElec;
	vector<double> eVisGam;
	double eVisAnn = 2.* 0.511;

	for(int branchIdx=0; branchIdx<m_nBr; branchIdx++)
	{
		//eVisGam[branchIdx] = 0;
		eVisGam.push_back(0);
		for (int gamIdx = 0; gamIdx < m_nGam; gamIdx++)
		{
			//if(m_eTruGam[branchIdx][gamIdx]>100)
				//eVisGam[branchIdx] = 0.9;
			//else 
			if(m_eTruGam[branchIdx][gamIdx]==0) break;
			//if(m_eTruGam[branchIdx][gamIdx]<1000)
			//if(m_eTruGam[branchIdx][gamIdx]>0)
			//eVisGam[branchIdx] += m_eTruGam[branchIdx][gamIdx];
            eVisGam[branchIdx] += EvisGamma(to_string(m_eTruGamStr[branchIdx][gamIdx]))*m_eTruGam[branchIdx][gamIdx] * (1+m_gammaScale);
			//eVisGam[branchIdx] += GetEVisGamma(m_eTruGam[branchIdx][gamIdx]);
					//eVisGam[branchIdx] = 0.9;
		}

		//if(m_eTruAlp[branchIdx]>0) {
			////std::cout << " ALPEVIS = " << m_eTruAlp[branchIdx]*dybEnergyModel::AlphaNL(m_eTruAlp[branchIdx]) << std::endl;
			//eVisGam[branchIdx] += m_eTruAlp[branchIdx]*dybEnergyModel::AlphaNL(m_eTruAlp[branchIdx]/3.);
		//}
    }
    //if(m_name == "B12") eVisGam[1] = 4.47167;

	/// apply scintillator NL
	for(int i=0; i<m_nBins; i++)
	{
		eTru     = m_binCenter[i];
        eVisElec = eTru;
        eVisElec = eTru * (electronQuench::ScintillatorNL(eTru)+electronCerenkov::getCerenkovPE(eTru));
        if(is_positron) eVisElec += 2*EvisGamma("511")*0.511 * (1+m_gammaScale) ;
        //if(is_positron) eVisElec += 2*0.467167;
        //cout << eTru << " " << electronQuench::getBirk1() << " " << electronQuench::ScintillatorNL(eTru)<< " " << electronCerenkov::getkC() << " " << eVisElec << endl;
		for(int branchIdx=0; branchIdx<m_nBr; branchIdx++)
		{
			eVis    = eVisElec + eVisGam[branchIdx];
            if(m_name=="B12" and branchIdx==2) eVis += 0.5; // 3-alpha branch
			eVisBck = eVis     + eVisAnn;
			//newBinSig = int(eVisSig/m_binWidth);
			//newBinBck = int(eVisBck/m_binWidth);
			newBinLow    = int((eVis   -m_eMin)/m_binWidth);
			newBinLowBck = int((eVisBck-m_eMin)/m_binWidth);
			newBinHig    = int((eVis   -m_eMin)/m_binWidth)+1;
			newBinHigBck = int((eVisBck-m_eMin)/m_binWidth)+1;
			bias         = (eVis    -m_eMin - newBinLow   *m_binWidth)/m_binWidth;
			biasBck      = (eVisBck -m_eMin - newBinLowBck*m_binWidth)/m_binWidth;
            //cout << newBinLow << " " << newBinHig << " " << bias << endl;
			//if(newBinSig<0) newBinSig=0;
			//if(newBinBck<0) newBinBck=0;
			//if(newBinSig>=m_nBins) newBinSig=m_nBins-1;
			//if(newBinBck>=m_nBins) newBinBck=m_nBins-1;
			//std::cout << " newBinSig = " << newBinSig << std::endl;
			//if(i>6500&&i<7000)
				//std::cout << branchIdx << "/" << i << ": newBinHigBck = " << newBinHigBck << std::endl;
			if(newBinLow<m_nBins)    m_eVis   [newBinLow   ] += (1-bias   )*m_eTru   [branchIdx][i];
			if(newBinLowBck<m_nBins) m_eVisBck[newBinLowBck] += (1-biasBck)*m_eTruBck[branchIdx][i];
			if(newBinHig<m_nBins)    m_eVis   [newBinHig   ] +=    bias    *m_eTru   [branchIdx][i];
			if(newBinHigBck<m_nBins) m_eVisBck[newBinHigBck] +=    biasBck *m_eTruBck[branchIdx][i];
		}
	}

    //for(int j=0; j<m_nBins; j++) {
    //    cout << m_eVis[j] << endl;
    //}

}

double junoSpectrum::GetChi2(int nDoF)
{
    ApplyScintillatorNL ();
    Normalize           ();
    double chi2 = 0;
    int rebin = m_nBins / m_nBinsData;
    double binWidthData = m_binWidth * rebin;
    m_nData = 0;
    for(int i=0; i < m_nBinsData; i++) {
        if(i*binWidthData<m_fitMin or binWidthData*i>m_fitMax-0.1) continue;
        if( m_eDataErr[i]!=0 ) {
            chi2 += pow( (m_eData[i] - m_eTheo[i])/m_eDataErr[i], 2); 
            m_nData++;
        }
    }
    cout << m_name << "  chi2: " << chi2 << " with nData : " << m_nData << endl;
	if(nDoF>0) chi2 /= double(m_nData - nDoF);
	return chi2;
}



double junoSpectrum::EvisGamma(string eTru)
{ 
    TFile file(junoParameters::gammaPdf_File.c_str(), "read");
    string pdfName = "gamma"+eTru+"keV";
    TH1D* gGammaPdf = (TH1D*)file.Get(pdfName.c_str());
    if(!gGammaPdf) cout << "No Such Pdf : " << pdfName << endl;
    //cout << "Read " << pdfName << " P.d.f ..." << endl;

    for(int i=0; i<gGammaPdf->GetNbinsX(); i++)  {
        m_pdf_eTru[i] = gGammaPdf->GetBinCenter(i+1);
        m_pdf_prob[i] = gGammaPdf->GetBinContent(i+1);
        if (m_pdf_prob[i] > 0) m_max_eTru = i;    // max energy cut ... 
    }
    file.Close();
    
    double numerator = 0.; double denominator = 0.;
    for(int iBin=1;  iBin<m_max_eTru; iBin++) {
        double E1 = m_pdf_eTru[iBin-1];
        double E2 = m_pdf_eTru[iBin];

        double prob1 = m_pdf_prob[iBin-1];
        double prob2 = m_pdf_prob[iBin];

        double fNL1 = electronQuench::ScintillatorNL(E1) + electronCerenkov::getCerenkovPE(E1);
        double fNL2 = electronQuench::ScintillatorNL(E2) + electronCerenkov::getCerenkovPE(E2);

        numerator   += ( prob1*E1*fNL1 + prob2*E2*fNL2 ) * (E2-E1) /2.;
        denominator += (prob1*E1 + prob2*E2) * (E2-E1)/ 2.;
    } 

    if(denominator ==0) { cout << " >> Error Happens While CalculateGammaNL <<<" << endl; return 0;}
    return numerator/denominator;

}













