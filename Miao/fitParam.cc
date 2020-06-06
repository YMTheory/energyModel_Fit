#include "fitParam.hh"

int main()
{
    SetStyle();
    
    //junoB12Data* m_b12data = new junoB12Data();
    ////std::cout << m_b12data->EvisGamma("4440") * 4.44<< endl;
    //m_b12data->SetParameters();
    //m_b12data->LoadData("./data/electron/B12.root");
    //cout << "Chi2: " << m_b12data->GetChi2(1) << endl;
    //m_b12data->Plot();

    //electronNLExperiment* junoENLExp = new electronNLExperiment();
    //junoENLExp->UpdateDataElectronNL();
    //junoENLExp->UpdateTheoElectronNL();
    //std::cout << junoENLExp->GetChi2() << endl;
    //junoENLExp->Plot();

    //electronNLChiFunction* junoEChiFCN = new electronNLChiFunction();
    //cout << junoEChiFCN->GetChiSquare() << endl;
    //junoEChiFCN->Plot();
    //junoEChiFCN->DrawContour(0,1);

    //gammaNLExperiment* junoGNLExp = new gammaNLExperiment();
    //junoGNLExp->UpdateTheoGammaNL();
    //cout << junoGNLExp->GetChi2(0) << endl;
    //junoGNLExp->Plot();

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

    junoNLChiFunction* junoNLChiFcn = new junoNLChiFunction();
    junoNLChiFcn->LoadData();
    ////cout << junoNLChiFcn->GetChi2() << endl; 
    junoNLChiFcn->GetChiSquare();
    junoNLChiFcn->Plot();
    //junoNLChiFcn->GridSearch();
    //junoNLChiFcn->ScanContour();
    //junoNLChiFcn->Plot();

    //junoC11Data* m_c11data = new junoC11Data();
    //m_c11data->SetParameters();
    //m_c11data->LoadData("./data/electron/C11.root");
    //cout << "C11 chi2: " << m_c11data->GetChi2() << endl;
    //m_c11data->Plot();
    //cout << m_c11data->EvisGamma("511")*0.511 << endl;
    


    return 1.0;
}
