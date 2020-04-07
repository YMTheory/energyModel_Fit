#ifndef _ELECTRONCERENKOV_H
#define _ELECTRONCERENKOV_H

#include <vector>

using namespace std;

class electronCerenkov
{
    public:
        electronCerenkov();
        ~electronCerenkov();
    
    public:
        void read_Cerenkov();
        
        double getCerenkovPE(double E);

    private:
        vector<double> Etrue;
        vector<double> Cerenkov;

};

#endif
