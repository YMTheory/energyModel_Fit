#include "fitParam.hh"

int main()
{
    SetStyle();

    gammaFitter*  gfitter = new gammaFitter();
    gfitter->Initialize();
    gfitter->Minimization();

    return 1.0;
}
