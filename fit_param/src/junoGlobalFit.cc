#include "junoGlobalFit.hh"

junoGammaData *junoGlobalFit::m_gammaData  = 0;

int    junoGlobalFit::m_nParameter = 4;
int    junoGlobalFit::m_nFreeParameter = 4;
double junoGlobalFit::m_parameters[4] ;

double junoGlobalFit::m_chi2       = 0;
double junoGlobalFit::m_chi2Gamma  = 0;
double junoGlobalFit::m_chi2Min    = -1;
double junoGlobalFit::m_bestFit     [20];
double junoGlobalFit::m_bestFitError[20];
double junoGlobalFit::m_covMatrix   [20][20];

junoGlobalFit::junoGlobalFit()
{
	m_gammaData  = new junoGammaData();
	cout<<"done1"<<endl;
}

junoGlobalFit::~junoGlobalFit()
{
    delete m_gammaData;
}

void junoGlobalFit::LoadData()  {
    std::cout << " >>> Start Loading Data >>> " << endl;
    m_gammaData->LoadData(junoParameters::gammaData_file);

}

double junoGlobalFit::GetChi2 (double maxChi2) {
    m_chi2 = 0;
    m_chi2Gamma = 0;

    if(junoParameters::fitGamma)  {
        m_chi2Gamma = m_gammaData->GetChi2();
        m_chi2     += m_chi2Gamma;
        if(maxChi2>0 && m_chi2>maxChi2) return 10000;
    }
    //m_chi2 += pow((junoGammaPeak::)) 

    return m_chi2;
}

void junoGlobalFit::SetParameters ()  {
    //cout << " >>> parameters Setting <<< " <<endl;
    junoEnergyModel::SetScintP0(m_parameters[0]);
    junoEnergyModel::SetScintP1(m_parameters[1]);
    junoEnergyModel::SetScintP2(m_parameters[2]);

    //std::cout << "Scint p0: " << junoEnergyModel::s_p0 << std::endl;
    //std::cout << "Scint p1: " << junoEnergyModel::s_p1 << std::endl;
    //std::cout << "Scint p2: " << junoEnergyModel::s_p2 << std::endl;
}

void junoGlobalFit::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag){
    //for (int i=0; i>m_nParameter; i++) {
    //    m_parameters[i] = par[i];
    //}
    //SetParameters  ();
    //cout << "Parameter Scan ---> p0: " << par[0] << " ---> p1: " << par[1] << " ---> p2: " << par[2] << endl;
    //cout << "Global Chi2 : " << GetChi2() << endl; 
    junoEnergyModel::SetScintP0(par[0]);
    junoEnergyModel::SetScintP1(par[1]);
    junoEnergyModel::SetScintP2(par[2]);
    fval = GetChi2 ();
}



void junoGlobalFit::Fit ()  {
    m_nParameter = 4;
    m_minuit = new TMinuit (m_nParameter);
    m_minuit->SetPrintLevel(junoParameters::fitPrintLevel);
    m_minuit->SetFCN(ChisqFCN);

    double arglist[4];
    int ierrflag = 0;

    m_minuit->mnexcm("CLEAR", arglist, 0, ierrflag);
    
    /// Scintillation Parameter 
    string scintPar[3] = {"kA", "kB", "kC"};
    double parMin[3] = {0, 1e-4, 0};
    double parMax[3] = {2, 1e-2, 2};
    double stepSize[3] = {0.1, 1e-5, 0.1};

    cout << "start point ---> " << junoParameters::p0_start << endl;
	m_minuit->mnparm(0,scintPar[0].c_str(), junoParameters::p0_start, stepSize[0],  parMin[0],parMax[0],ierrflag);
	m_minuit->mnparm(1,scintPar[1].c_str(), junoParameters::p1_start, stepSize[1],  parMin[1],parMax[1],ierrflag);
	m_minuit->mnparm(2,scintPar[2].c_str(), junoParameters::p2_start, stepSize[2],  parMin[2],parMax[2],ierrflag);

    // pull terms for gamma relative energy scale
    m_minuit->mnparm(3, "gammaScale", junoParameters::gamScale, 0.1*junoParameters::gamScaleError, 0.5, 1.5, ierrflag);

    // fix parameters 
    if(junoParameters::fixScintP0)   m_minuit->FixParameter(0);
    if(junoParameters::fixScintP1)   m_minuit->FixParameter(1);
    if(junoParameters::fixScintP2)   m_minuit->FixParameter(2);
    if(junoParameters::fixGamScale)  m_minuit->FixParameter(3);

    m_nParameter = m_minuit->GetNumPars();

	// Minimization strategy
	// 1 standard;	2 try to improve minimum (slower)
	arglist[0] = 2;
	m_minuit->mnexcm("SET STR", arglist, 1, ierrflag);

    arglist[0] = 200000;     /// Max Calls
    arglist[1] = 0.01;      // tolerance

    m_minuit->mnexcm("MIGrad", arglist, 2, ierrflag);

	double min, edm, errdef;
	int nvpar, nparx, icstat;
	m_minuit->mnstat(m_chi2Min, edm, errdef, nvpar, nparx, icstat);

	for(int i=0; i<m_nParameter; i++)
	{
		m_minuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
	}
	for (int i = 0; i < m_nParameter; i++)
		m_parameters[i] = m_bestFit[i];
	SetParameters();  
	
	m_nFreeParameter = m_minuit->GetNumFreePars() - 1;
	m_minuit->mnemat(&m_covMatrix[0][0],20);
	delete m_minuit;
	std::cout << " =============================== " << std::endl;
	std::cout << " chi2 minimum: " << m_chi2Min << std::endl;
	std::cout << " =============================== " << std::endl;
}
