#ifndef _PRIMARYELECTRONDISTRIBUTION_H
#define _PRIMARYELECTRONDISTRIBUTION_H

#include <vector>

using namespace std;

class primaryElectronDistribution
{
    public:
        primaryElectronDistribution();
        ~primaryElectronDistribution();

    public:
        void read_distribution(string source);

        int getSize() { return Etrue.size(); }
        double getEtrue(int num);
        int getCount(int num);

        
    private:
        vector<double> Etrue;
        vector<double> prm_count;

};

#endif