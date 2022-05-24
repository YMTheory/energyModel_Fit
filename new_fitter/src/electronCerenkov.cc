#include "electronCerenkov.hh"
#include "junoParameters.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <TMath.h>
#include <TGraph.h>
#include <TFile.h>

using namespace std;

double electronCerenkov::m_kC = 3;
//double electronCerenkov::m_energyScale = junoParameters::m_energyscale;
double electronCerenkov::m_energyScale = 1409.84; //3134.078/2.223;
bool electronCerenkov::m_LoadCerenkov = false;
double electronCerenkov::m_A1 = 0;
double electronCerenkov::m_A2 = -7.34716;  //-7.90041506;
double electronCerenkov::m_A3 = 15.5519;   //13.84029077;
double electronCerenkov::m_A4 = 0.0267155; //0.03641905;
double electronCerenkov::m_E0 = 0.2;

//double electronCerenkov::m_p0 = -0.410;
//double electronCerenkov::m_p1 =  0.42;
//double electronCerenkov::m_p2 = -0.0165;
//double electronCerenkov::m_p3 = 3.133;
//double electronCerenkov::m_p4 = 0.0225;
double electronCerenkov::m_p0 = 0.4;
double electronCerenkov::m_p1 = 0.122;    
double electronCerenkov::m_p2 = 0.573;      
double electronCerenkov::m_p3 = 0.816455 * 8.701754386 * 20;    
double electronCerenkov::m_p4 =  0.178972;      
//double electronCerenkov::m_p0 = 0.4;
//double electronCerenkov::m_p1 = 0.122;    
//double electronCerenkov::m_p2 = 0.573;      
//double electronCerenkov::m_p3 = 0.816455 * 8.701754386;    
//double electronCerenkov::m_p4 =  0.178972;      

vector<double> electronCerenkov::m_Etrue;
vector<double> electronCerenkov::m_Cerenkov;

vector<double> electronCerenkov::m_Etrue_local;
vector<double> electronCerenkov::m_Cerenkov_local;

double electronCerenkov::m_E[m_nData];
double electronCerenkov::m_nonl[m_nData];

TGraph* electronCerenkov::gNPE_elec = new TGraph();
TGraph* electronCerenkov::gNPE_elec_local = new TGraph();

electronCerenkov::electronCerenkov()
{}

electronCerenkov::~electronCerenkov()
{;}

void electronCerenkov::LoadCerenkov()
{
    cout << " >>> Loading Electron Cerenkov Shape <<< " << endl;
    ifstream in;
    in.open(junoParameters::cerenkovNL_File.c_str());
    if(!in){
        cout << " >>> Fail to Open Cerenkov File !! <<< " << endl;
    }
    string line;

    double tmp_Edep, tmp_rel, tmp_abs; int index = 0;
    while(getline(in,line)){
        istringstream ss(line);
        //ss >> tmp_Edep >> tmp_rel >> tmp_abs ;
        ss >> tmp_Edep  >> tmp_abs ;
        m_Etrue.push_back(tmp_Edep);
        //m_Etrue.push_back(tmp_Edep/1000.);
        m_Cerenkov.push_back(tmp_abs);
        gNPE_elec->SetPoint(index, tmp_Edep, tmp_abs);
        index++;
    }

    in.close();

    TFile* finCer = new TFile("/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/pyfitter/data/Ncer_local.root", "read");
    TGraph* gNcer = (TGraph*)finCer->Get("Ncer");
    for(int i=0; i<gNcer->GetN(); i++) {
        m_Etrue_local.push_back(gNcer->GetPointX(i));
        m_Cerenkov_local.push_back(gNcer->GetPointY(i));
        gNPE_elec_local->SetPoint(i, gNcer->GetPointX(i), gNcer->GetPointY(i));
    }

    delete gNcer;
    delete finCer;

    m_LoadCerenkov = true;
}


double electronCerenkov::getCerenkovPE(double E)
{
    if(!m_LoadCerenkov)   LoadCerenkov();

    if(m_Cerenkov.size() == 0) {
        cout << " >>> No Data in Cerenkov Vector <<< " << endl; return -1;
    } else if (m_Cerenkov.size() != m_Etrue.size()){
        cout << " >>> Cerenkov Vector Length are Different !! <<< " << endl; 
    } else {

        // get Cerenkov PE
        int num = m_Cerenkov.size();
        for(int i=1; i<num; i++){
            if(m_Etrue[i-1]<=E and m_Etrue[i]>=E){  
                return m_kC*((m_Cerenkov[i]*(E-m_Etrue[i-1])+m_Cerenkov[i-1]*(m_Etrue[i]-E))/(m_Etrue[i]-m_Etrue[i-1]))/m_energyScale/E;
                //return m_kC * m_Cerenkov[i-1]/m_energyScale/E; 
            }
        }
        cout << E <<  " >>> Energy Beyond Range !! <<< " << endl; return -1;
    }
}


double electronCerenkov::getSimCerPE(double E)
{
    if(!m_LoadCerenkov)   LoadCerenkov();

    //if(m_Cerenkov.size() == 0) {
    //    cout << " >>> No Data in Cerenkov Vector <<< " << endl; return -1;
    //} else if (m_Cerenkov.size() != m_Etrue.size()){
    //    cout << " >>> Cerenkov Vector Length are Different !! <<< " << endl; 
    //} else {

    //    // get Cerenkov PE
    //    int num = m_Cerenkov.size();
    //    for(int i=1; i<num; i++){
    //        if(m_Etrue[i-1]<=E and m_Etrue[i]>=E){  
    //            return m_kC*((m_Cerenkov[i]*(E-m_Etrue[i-1])+m_Cerenkov[i-1]*(m_Etrue[i]-E))/(m_Etrue[i]-m_Etrue[i-1]));
    //            //return m_kC * m_Cerenkov[i-1]/m_energyScale/E; 
    //        }
    //    }
    //    cout << E <<  " >>> Energy Beyond Range !! <<< " << endl; return -1;
    //}

    //return m_kC*gNPE_elec->Eval(E, 0, "S");
    return m_kC*gNPE_elec->Eval(E);
}

double electronCerenkov::getLocalSimCerPE(double E)
{
    if (!m_LoadCerenkov) LoadCerenkov();
    return m_kC * gNPE_elec_local->Eval(E);

}


void electronCerenkov::Plot()
{
    cout << " >>> Plot Cerenkov Curve <<< " << endl;
    for (int i=0; i<m_nData; i++) {
        m_E[i] = m_Etrue[i];
        m_nonl[i] = getCerenkovPE(m_E[i]);
    }

    TFile* file = new TFile(junoParameters::cerenkov_outFile.c_str(), "recreate");
    TGraph* gCerenkovNL = new TGraph(1000, m_E, m_nonl);
    gCerenkovNL->SetLineColor(kBlue+1);
    gCerenkovNL->SetMarkerColor(kBlue+1);
    gCerenkovNL->SetMarkerSize(0.2);
    gCerenkovNL->Write();
    file->Close();
}

double electronCerenkov::getAnaCerPE(double E)
{
    double x = TMath::Log(1+E/m_E0);
    double npe = (m_A2*x*x + m_A3*x*x*x) * (1/E + m_A4) * E   ;
    return npe;

}



double electronCerenkov::getNewAnaCerPE(double E)
{
    if (E < m_E0)
        return 0;
    else{
        //E = E - 0.2;   // 0.2 MeV threshold
        // E0 should be a free parameter?
        E = E - m_E0;
        double NC = (m_p3*E*E) / (m_p4+m_p0*E+m_p1*TMath::Exp(-m_p2*E));
        return NC;
    }
}


double electronCerenkov::getNewAnaCerPE1(double E)
{
    if (E < m_E0)
        return 0;
    else{
        //E = E - 0.2;   // 0.2 MeV threshold
        // E0 should be a free parameter?
        E = E - m_E0;
        double NC = (m_p0*E*E) / (E+m_p1*TMath::Exp(-m_p2*E));
        return NC;
    }
}



double electronCerenkov::getCerPE(double E) {
    if (junoParameters::cerenkovMode == "kSimulationCer" )
        return getSimCerPE(E);
    else if (junoParameters::cerenkovMode == "kSimulationCerLocal" ) {
        return getLocalSimCerPE(E);
    } else if (junoParameters::cerenkovMode == "kAnalyticalCer" ) {
        return getAnaCerPE(E);
    } else if (junoParameters::cerenkovMode == "kAnalyticalNewCer") {
        return getNewAnaCerPE(E);
    } else if (junoParameters::cerenkovMode == "kAnalyticalNewCer1") {
        return getNewAnaCerPE1(E);
    }
}


