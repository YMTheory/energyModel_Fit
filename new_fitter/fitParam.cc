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

    //gammaData* cs137 = new gammaData("nFe56", 700, 1100, 100);
    //cs137->LoadData();
    //cs137->calcGammaResponse();
    //cout << cs137->GetChi2();
    
    //gammaResponse* cs137 = new gammaResponse("Cs137", 2000, 0, 15000);
    //cs137->LoadData();
    //cs137->calcGamResponse();
    //gammaResponse* Mn54 = new gammaResponse("Mn54", 2000, 0, 15000);
    //Mn54->LoadData();
    //Mn54->calcGamResponse();
    //gammaResponse* K40 = new gammaResponse("K40", 2000, 0, 15000);
    //K40->LoadData();
    //K40->calcGamResponse();
    //gammaResponse* nH = new gammaResponse("nH", 2000, 0, 15000);
    //nH->LoadData();
    //nH->calcGamResponse();
    //gammaResponse* Tl208 = new gammaResponse("Tl208", 2000, 0, 15000);
    //Tl208->LoadData();
    //Tl208->calcGamResponse();
    //gammaResponse* nC12 = new gammaResponse("nC12", 2000, 0, 15000);
    //nC12->LoadData();
    //nC12->calcGamResponse();
    //gammaResponse* O16 = new gammaResponse("O16", 2000, 0, 15000);
    //O16->LoadData();
    //O16->calcGamResponse();
    //gammaResponse* nFe56 = new gammaResponse("nFe56", 2000, 0, 15000);
    //nFe56->LoadData();
    //nFe56->calcGamResponse();

    //cout << cs137->GetName() << " " << cs137->GetNonlData() << " " << cs137->GetNonlCalc() << endl;
    //cout << Mn54->GetName() << " " << Mn54->GetNonlData() << " " << Mn54->GetNonlCalc() << endl;
    //cout << K40->GetName() << " " << K40->GetNonlData() << " " << K40->GetNonlCalc() << endl;
    //cout << nH->GetName() << " " << nH->GetNonlData() << " " << nH->GetNonlCalc() << endl;
    //cout << Tl208->GetName() << " " << Tl208->GetNonlData() << " " << Tl208->GetNonlCalc() << endl;
    //cout << nC12->GetName() << " " << nC12->GetNonlData() << " " << nC12->GetNonlCalc() << endl;
    //cout << O16->GetName() << " " << O16->GetNonlData() << " " << O16->GetNonlCalc() << endl;
    //cout << nFe56->GetName() << " " << nFe56->GetNonlData() << " " << nFe56->GetNonlCalc() << endl;

    return 1.0;
}
