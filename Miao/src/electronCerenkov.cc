#include "electronCerenkov.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <TMath.h>
#include <TGraph.h>

using namespace std;

electronCerenkov::electronCerenkov()
{;}

electronCerenkov::~electronCerenkov()
{;}

void electronCerenkov::read_Cerenkov()
{
    ifstream in;
    in.open("./data/absolutePE.txt");
    if(!in){
        cout << " >>> Fail to Open Cerenkov File !! <<< " << endl;
    }
    string line;

    double tmp_Edep, tmp_rel, tmp_abs;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_Edep >> tmp_rel >> tmp_abs ;
        Etrue.push_back(tmp_Edep/1000.);
        Cerenkov.push_back(tmp_abs);
    }

    in.close();
}


double electronCerenkov::getCerenkovPE(double E)
{
    read_Cerenkov();

    if(Cerenkov.size() == 0) {
        cout << " >>> No Data in Cerenkov Vector <<< " << endl; return -1;
    } else if (Cerenkov.size() != Etrue.size()){
        cout << " >>> Vector Length are Different !! <<< " << endl; 
    } else {

        // get Cerenkov PE
        int num = Cerenkov.size();
        for(int i=1; i<num; i++){
            if(Etrue[i-1]<=E and Etrue[i]>E){  return Cerenkov[i]; }
        }
        cout << E <<  " >>> Energy Beyond Range !! <<< " << endl; return -1;
    }
}
