#include "fitParam.hh"

int main()
{
    SetStyle();

    //gammaResol * gResol = new gammaResol("Co60", 100, 700, 1100);
    //gResol->check_nonl() ;
    //gResol->calcGammaNPE();
    //cout << "initial chi2 value: " << gResol->GetChi2() << endl;
    //gResol->Plot();
    //delete gResol;
    //cout << "Predicted Resolution : " << gResol->GetResolPred() << endl; 


    junoNLChiFunction* junoNLChi2 = new junoNLChiFunction();
    junoNLChi2->LoadData();
    junoNLChi2->GetChiSquare();
    junoNLChi2->Plot();

    return 1.0;
}
