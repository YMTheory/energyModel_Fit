#include "junoEnergyModel.hh"

double junoEnergyModel::s_p0   = 0.98;
double junoEnergyModel::s_p1   = 6.5e-3;
double junoEnergyModel::s_p2   = 1.;
double junoEnergyModel::s_p3   = 0.;

double junoEnergyModel::s_kB    = 6.5e-3;
double junoEnergyModel::s_cer     = 1.;
double junoEnergyModel::s_rad     = 0.;

bool   junoEnergyModel::s_isLoaded = false;
double junoEnergyModel::s_kBResid  = 0;

double junoEnergyModel::s_energySamples         [s_nSamples] = {0};
double junoEnergyModel::s_quenchingShape1[s_nKb][s_nSamples] = {0};
double junoEnergyModel::s_cerenkovShape         [s_nSamples] = {0};

double* junoEnergyModel::s_quenchingShape1_lowKb = &s_quenchingShape1[0][0];
double* junoEnergyModel::s_quenchingShape1_higKb = &s_quenchingShape1[0][0];

vector<double> junoEnergyModel::m_energySamples;


junoEnergyModel::junoEnergyModel()
{;}

junoEnergyModel::~junoEnergyModel()
{;}

void junoEnergyModel::Load(){
	std::cout << " ----> Loading quenching+cherenkov shapes " << std::endl;
    // load quenching shapes ...
    TFile* quenchingFile = new TFile("data/Quench.root", "read");
    for(int kbIdx=10; kbIdx<s_nKb; kbIdx++)  {
        //cout << kbIdx << endl;
        stringstream ss; ss << kbIdx;
        TString name1 = "kB"+ss.str();
        TGraph* quench1G = (TGraph*)quenchingFile->Get(name1);
        double* quench1 = quench1G->GetY();

        for(int sampleIdx=0; sampleIdx<s_nSamples; sampleIdx++)
        {
			s_quenchingShape1[kbIdx][sampleIdx] = quench1[sampleIdx];
        }
        delete quench1G;
    }
	quenchingFile->Close();
    delete quenchingFile;

    // load Cerenkov shapes ...
    ifstream infileC("./data/CerenkovPE.txt");
    TGraph cerenkovG(0);
    double xC, rC, aC;
    for(Int_t i=0; i!=799; i++){
        infileC >> xC >> rC >> aC;
        cerenkovG.SetPoint(i, xC/1000., aC);
    }
    infileC.close();
    for(int sampleIdx=0; sampleIdx<s_nSamples; sampleIdx++)  {
        double energy = sampleIdx * s_samplingResol;
        s_cerenkovShape[sampleIdx] = cerenkovG.Eval(energy);
    }
	s_isLoaded = true;
}

void junoEnergyModel::Update(){
    if(!s_isLoaded) Load();
    if(s_kB==0) return;
    int kBIdx = int(s_kB*1e4); //cout << "kBidx" << kBIdx << endl; 
	s_kBResid = kBIdx+1 - s_kB*1e4; 
	s_quenchingShape1_lowKb = &s_quenchingShape1[kBIdx]  [0];
	s_quenchingShape1_higKb = &s_quenchingShape1[kBIdx+1][0];
}

double junoEnergyModel::ScintillatorNL(double eTru) {

    //double norm = 1.;
    //if(s_normEnergy>0){
    //    norm = ScintillatorShape(s_normEnergy);
    //}
    double nl = ScintillatorShape(eTru);
    return nl;
}

double junoEnergyModel::ScintillatorShape(double eTru)  {
    if(!s_isLoaded) Load();
    Update();
    return PhysicsScintillator(eTru);
}

double junoEnergyModel::PhysicsScintillator(double eTru)  {
    if(s_kB==0) return 1.0;

    int idx = int(eTru/s_samplingResol+0.5);

    double cerenkNL = s_cer * s_cerenkovShape[idx] / 1481.06 / eTru;

    double quenchNL  =  s_p0* ( s_kBResid    *s_quenchingShape1_lowKb[idx] 
                    +(1-s_kBResid) *s_quenchingShape1_higKb[idx] );
  
    // check parameters ... 
    //std::cout << " -----------> Parameters <---------" << std::endl;
    //std::cout << " -----> idx : "      << idx      << std::endl;
    //std::cout << " -----> s_cer : "    << s_cer    << std::endl;
    //std::cout << " -----> s_p0 : "     << s_p0     << std::endl;
    //std::cout << " -----> cerenkNL : " << cerenkNL << std::endl;
    //std::cout << " -----> scintNL : "  << quenchNL << std::endl;

    return quenchNL+cerenkNL;
}

double junoEnergyModel::QuenchNL (double eTru) {
    
    int idx = int(eTru/s_samplingResol+0.5);
    double quenchNL  =  s_p0* ( s_kBResid * s_quenchingShape1_lowKb[idx] 
                    +(1-s_kBResid) *s_quenchingShape1_higKb[idx] );
    return quenchNL;
}

double junoEnergyModel::CerenkovNL (double eTru) {
    int idx = int(eTru/s_samplingResol+0.5);
    double cerenkNL = s_cer * s_cerenkovShape[idx] / 1481.06 / eTru;
    return cerenkNL; 
}

void junoEnergyModel::SaveCurves() {
    cout << " kA = " << junoEnergyModel::s_p0    << endl;
    cout << " kB    = " << junoEnergyModel::s_p1    << endl;
    cout << " kC    = " << junoEnergyModel::s_p2    << endl;
    
    TGraph eleS = DrawElectronScintNL();
    TGraph eleQ = DrawElectronQuenchNL();
    TGraph eleC = DrawElectronCerenkovNL();
    //TGraph gamS = DrawGammaScintNL   ();
    string name = "output/curves/curves_nonl.root";
    TFile* nlFile = new TFile(name.c_str(),"recreate");
    eleS.Write("electronScintNL");
    eleQ.Write("electronQuenchNL");
    eleC.Write("electronCerenkovNL");
    nlFile->Close();
    delete nlFile;
}

TGraph junoEnergyModel::DrawElectronScintNL(int nSamples, double eMax)
{
    TGraph gr(0);
    gr.SetLineColor(kBlue-2);
    m_energySamples.clear();
    vector<double> nl = SampleElectronScintNL(nSamples,eMax);
    for(int i=0; i<m_energySamples.size(); i++)
    {
        gr.SetPoint(i,m_energySamples[i],nl[i]);
    }
    return gr;

}

TGraph junoEnergyModel::DrawElectronQuenchNL(int nSamples, double eMax) 
{
    TGraph gr(0);
    gr.SetLineColor(kBlue-2);
    m_energySamples.clear();
    vector<double> nl = SampleElectronQuenchNL(nSamples,eMax);
    for(int i=0; i<m_energySamples.size(); i++)
    {
        gr.SetPoint(i,m_energySamples[i],nl[i]);
    }
    return gr;
}


TGraph junoEnergyModel::DrawElectronCerenkovNL(int nSamples, double eMax) 
{
    TGraph gr(0);
    gr.SetLineColor(kBlue-2);
    m_energySamples.clear();
    vector<double> nl = SampleElectronCerenkovNL(nSamples,eMax);
    for(int i=0; i<m_energySamples.size(); i++)
    {
        gr.SetPoint(i,m_energySamples[i],nl[i]);
    }
    return gr;
}


vector<double> junoEnergyModel::SampleElectronScintNL(int nSamples, double eMax)  {
    
    m_energySamples.clear();
    vector<double> electronScintNL;
    double deltaE = eMax/nSamples;
    double eTru;
    for(int i=0; i<nSamples; i++){
        eTru = (i+1) * deltaE;
        if(eTru > eMax)  break;
        m_energySamples.push_back(eTru);
        electronScintNL.push_back(ScintillatorNL(eTru));
    }
    return electronScintNL;
}


vector<double> junoEnergyModel::SampleElectronQuenchNL (int nSamples, double eMax) {
    m_energySamples.clear();
    vector<double> electronQuenchNL;
    double deltaE = eMax/nSamples;
    double eTru;
    for(int i=0; i<nSamples; i++){
        eTru = (i+1) * deltaE;
        if(eTru > eMax)  break;
        m_energySamples.push_back(eTru);
        electronQuenchNL.push_back(QuenchNL(eTru));
    }
    return electronQuenchNL;
}


vector<double> junoEnergyModel::SampleElectronCerenkovNL (int nSamples, double eMax) {
    m_energySamples.clear();
    vector<double> electronCerenkovNL;
    double deltaE = eMax/nSamples;
    double eTru;
    for(int i=0; i<nSamples; i++){
        eTru = (i+1) * deltaE;
        if(eTru > eMax)  break;
        m_energySamples.push_back(eTru);
        electronCerenkovNL.push_back(CerenkovNL(eTru));
    }
    return electronCerenkovNL;
}






