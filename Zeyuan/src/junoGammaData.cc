#include "junoGammaData.hh"

junoGammaData::junoGammaData()
{
    cout << " >>>>>> Initializing gamma peaks " << endl;
}

void junoGammaData::AddPeak (string name, string pdfName, double eTruSingle, double eTruTotal)
{
    std::cout << "  adding " << name << " peak " << std::endl;
    junoGammaPeak peak (name.c_str(), pdfName.c_str(), eTruSingle, eTruTotal);
    m_data.push_back (peak) ;

}

void junoGammaData::LoadData(string fileName)  {
    m_data.clear();
	///     --name--        --pdf-- --true E----
	///                                     total single
    AddPeak("^{137}Cs"     ,"gammaCs137",  0.6617, 0.6617);
    AddPeak("^{54}Mn"      ,"gammaMn54",   0.8348, 0.8348);
    AddPeak("^{68}Ge"      ,"gammaGe68",   1.0220, 1.0220);
    AddPeak("^{40}K"       ,"gammaK40",    1.4608, 1.4608);
    AddPeak("n-H"          ,"gammanH",     2.2233, 2.2233);
    //AddPeak("n-^{12}C"     ,"gammanC12",   4.9450, 4.9450);
    AddPeak("^{60}Co"      ,"gammaCo60",   2.5060, 2.5060);
    
    std::cout  << "Reading gamma data from " << fileName << std::endl;
	ifstream infile(fileName.c_str());
    if(!infile) { std::cout << "Fail to Open " << fileName << std::endl; }
    double eVis, eVisError;
    int i=0;
    while (infile >> eVis >> eVisError)  {
        if(i>=m_data.size() )
			std::cout << " ERROR: Data file contains more gamma peaks than initialized" << std::endl;
        m_data[i].SetEVis       (eVis);
        m_data[i].SetEVisError  (eVisError);  
        i++;
    }
	if(i!=m_data.size())
		std::cout << " ERROR: Data file contains less gamma peaks than initialized" << std::endl;
}

double junoGammaData::GetChi2 ( int nDoF )  {
    double chi2 = 0;
    vector<junoGammaPeak>::iterator peakItr = m_data.begin();
    for(; peakItr!=m_data.end(); peakItr++)
    {
        peakItr->UpdateDataNL();
        //peakItr->UpdateTheoNL();
        chi2 += peakItr->GetChi2();
    }
    if (nDoF>0)
    {
        chi2 /= double( m_data.size() - nDoF);
        m_nData = m_data.size();
    }
    return chi2;
}


void junoGammaData::GenToyMC()
{
	//std::cout << " ------> Generating " << junoParameters::nToy << " gamma toy MC samples " << std::endl;
    

}

TGraphErrors junoGammaData::PlotPeaks (int type)
{
    cout << " ---> Initializing Gamma Data ... " << endl;
    TGraphErrors gr(0);
	gr.SetLineColor  (kBlue+2);
	gr.SetMarkerColor(kBlue+2);
	gr.SetMarkerStyle(10);
	gr.SetMarkerSize (0.6);
    vector<junoGammaPeak>::iterator peakItr = m_data.begin();
    int pointIdx = 0;
    for(;peakItr!=m_data.end(); peakItr++)
    {
		std::cout << peakItr->GetName() << std::endl;
		double energy   = peakItr->GetETruSingle();
        double error    = peakItr->GetEVisError();
        double nl       = 1;
        peakItr->UpdateDataNL();
        if(type==1) nl = peakItr->GetDataScintNL();
        if(type==2) nl = peakItr->GetTheoScintNL();
		std::cout << " nl = " << nl << std::endl;
		std::cout << " error = " << error << std::endl;
		gr.SetPoint     (pointIdx,energy,nl);
		gr.SetPointError(pointIdx,0,nl*error);
		pointIdx++;
        
    }

    return gr;
}

void junoGammaData::Plot (bool writeToFile)
{
	cout << " ---> Intializing Gamma data " << endl;
    TGraphErrors* peaksScint    = new TGraphErrors;
	TGraphErrors* peaksScintR  = new TGraphErrors;
    peaksScint    ->SetLineColor(kBlue+1);
	peaksScint    ->SetMarkerColor(kBlue+1);
	peaksScint    ->SetMarkerStyle(20);
	peaksScint    ->SetMarkerSize(1.1);
	peaksScintR ->SetLineColor  (kBlue+1);
	peaksScintR ->SetMarkerColor(kBlue+1);
	peaksScintR ->SetMarkerSize(1.1);
    vector<junoGammaPeak>::iterator peakItr = m_data.begin();
	int pointIdx    = 0;
	for(;peakItr!=m_data.end();peakItr++)
	{
		peakItr->UpdateTheoNL();
		peakItr->UpdateDataNL();
		double eVis     = peakItr->GetEVis();
		double eTru     = peakItr->GetETruTotal();
		double error    = peakItr->GetEVisError();
        double nlScint  = peakItr->GetDataScintNL();
		double eVisPred = peakItr->GetTheoScintNL()*eTru;
        peaksScint->SetPoint(pointIdx, eTru, nlScint);
		peaksScint ->SetPointError(pointIdx,0,nlScint*error);
		peaksScintR->SetPoint     (pointIdx,eTru,eVis/eVisPred);
		peaksScintR->SetPointError(pointIdx,0,error);

    }

	string outname = "output/gammas_peaks.root";
	TFile* gammaFile = new TFile(outname.c_str(),"recreate");
	peaksScint  ->Write("peaksSingle")   ;
	peaksScintR ->Write("peaksScintR")  ;
	gammaFile->Close();
	delete gammaFile;
}

