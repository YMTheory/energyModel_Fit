#include "fitParam.hh"

int main()
{
    SetStyle();

    // gamma calibration source
    //gammaFitter*  gfitter = new gammaFitter();
    //gfitter->Initialize();
    //gfitter->Minimization();
    //gfitter->Plot();
    
    //junoSpectrum* B12 = new junoSpectrum(14400, 100,
    //                                     3, 1,
    //                                     0, 15,
    //                                     4, 14,
    //                                     "B12");
    //B12->LoadData();
    //B12->GetChi2();
    //B12->Plot();

    // beta spectum
    //spectrumFitter* fitter = new spectrumFitter();
    //fitter->Initialize();
    //fitter->Minimization();

    junoNLChiFunction* fitter = new junoNLChiFunction();
    fitter->LoadData();
    fitter->GetChiSquare();
    fitter->Plot();


    return 1.0;
}
