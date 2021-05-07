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

double electronCerenkov::m_kC = 1;
//double electronCerenkov::m_energyScale = junoParameters::m_energyscale;
double electronCerenkov::m_energyScale = 1409.84; //3134.078/2.223;
bool electronCerenkov::m_LoadCerenkov = false;
double electronCerenkov::m_A1 = 0;
double electronCerenkov::m_A2 = -7.90041506;
double electronCerenkov::m_A3 = 13.84029077;
double electronCerenkov::m_A4 = 0.03641905;
double electronCerenkov::m_E0 = 0.2;

vector<double> electronCerenkov::m_Etrue;
vector<double> electronCerenkov::m_Cerenkov;

double electronCerenkov::m_E[m_nData];
double electronCerenkov::m_nonl[m_nData];

TGraph* electronCerenkov::gNPE_elec = new TGraph();

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
    return (m_A1*x + m_A2*x*x + m_A3*x*x*x) * (1/E + m_A4) * E   ;

}


double electronCerenkov::getCerPE(double E) {
    if (junoParameters::cerenkovParameterization == kSimulationCer)
        return getSimCerPE(E);
    else if (junoParameters::cerenkovParameterization == kAnalyticalCer)
        return getAnaCerPE(E);
}


