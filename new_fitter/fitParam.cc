#include "fitParam.hh"

int main()
{
    SetStyle();

    gammaFitter*  gfitter = new gammaFitter();
    gfitter->Initialize();
    gfitter->Minimization();
    gfitter->Plot();

    return 1.0;
}
