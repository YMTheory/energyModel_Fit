#include "primaryElectronDistribution.hh"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

primaryElectronDistribution::primaryElectronDistribution()
{;}

primaryElectronDistribution::~primaryElectronDistribution()
{;}

void primaryElectronDistribution::read_distribution(string source)
{
    Etrue.clear();
    prm_count.clear();

    string name1 = "./data/naked_gamma/primary_";
    string name2 = ".txt";
    string filename = name1 + source + name2;

    ifstream in;
    in.open(filename.c_str());
    if(!in){ cout << "Error Open Primary e+- File" << endl; }
    string line;
    
    double tmp_E, tmp_count;
    while(getline(in,line)){
        istringstream ss(line);
        ss >> tmp_E >> tmp_count;
        Etrue.push_back(tmp_E);
        prm_count.push_back(tmp_count);
    }
    in.close();

}

double primaryElectronDistribution::getEtrue(int num)
{
    if(num > Etrue.size()) { cout << " >>> Error in Primary e+- Dist! Beyong File Range !! <<<" << endl;  return -1;}
    else
        return Etrue[num] ;
}

int primaryElectronDistribution::getCount(int num)
{
    if(num > prm_count.size()) { cout << " >>> Error in Primary e+- Dist! Beyong File Range !! <<<" << endl;  return -1;}
    else
        return prm_count[num];
    
}
