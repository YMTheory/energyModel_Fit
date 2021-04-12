#include "fitParam.hh"

int main()
{
    SetStyle();

    // electron fake nonl fitting :
    //electronFitter* efitter = new electronFitter();
    //efitter->Initialize();
    //efitter->Minimization();
    //efitter->Plot();

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
    //

    //electronResponse::FitPlot();

    //junoSpectrum* junoB12data = new junoSpectrum(14400, 100, 3, 2,
    //                         0, 15, 0, 15, "histogram", "B12");
    //junoB12data->LoadData();
    //cout << junoB12data->GetChi2() <<endl;
    //junoB12data->Plot();

    //junoB12* b12 = new junoB12();
    //b12->Initialize();
    //b12->LoadDataSpec();
    //b12->LoadTheoSpec();
    //cout << b12->GetChi2() << endl;
    //b12->Plot();

    return 1.0;
}
