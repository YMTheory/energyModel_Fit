#include "fitParam.hh"

int main()
{
    SetStyle();
    
    //electronNLExperiment* junoENLExp = new electronNLExperiment();
    //junoENLExp->UpdateDataElectronNL();
    //junoENLExp->UpdateTheoElectronNL();
    //std::cout << junoENLExp->GetChi2() << endl;
    //junoENLExp->Plot();

    BetaPrediction* pred = new BetaPrediction();
    pred->Plot();

    //electronNLChiFunction* junoEChiFCN = new electronNLChiFunction();
    //cout << junoEChiFCN->GetChiSquare() << endl;
    //junoEChiFCN->Plot();
    //junoEChiFCN->DrawContour(0,1);

    //gammaNLExperiment* junoGNLExp = new gammaNLExperiment();
    //junoGNLExp->UpdateTheoGammaNL();
    //cout << junoGNLExp->GetChi2(0) << endl;

    //gammaNLChiFunction* junoGChiFCN = new gammaNLChiFunction();
    //junoGChiFCN->GetChiSquare();
    //junoGChiFCN->Plot();
    //junoGChiFCN->DrawContour(0,1);

    //electronResExperiment* junoElecResExp = new electronResExperiment();
    //cout << "chi2: " << junoElecResExp->GetChi2() << endl;
    //junoElecResExp->PlotC11();

    //electronResChiFunction* junoEResChiFCN = new electronResChiFunction();
    //junoEResChiFCN->GetChiSquare();
    //junoEResChiFCN->Plot();

    return 1.0;
}
