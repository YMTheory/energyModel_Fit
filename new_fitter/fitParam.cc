#include "fitParam.hh"

int main()
{
    //SetStyle();

    
    //junoChiFunction* fitter = new junoChiFunction();
    //fitter->LoadData();
    //fitter->GetChiSquare();
    /////////fitter->setFixParId(1);
    /////////for (int i=0; i<100; i++) {
    /////////    double kB = ( 6.31677e-03 - 6.44315e-4 ) + 6.44315e-4 / 50 * i;
    /////////    fitter->setFixParValue(kB);
    /////////    fitter->GetChiSquare();
    /////////}
    //fitter->Plot();


    //junoB12_simplified* b12 = new junoB12_simplified(100, 4000, 18000);
    //b12->Initialize();
    ////b12->LoadDataSpec();
    //////b12->LoadTheoSpec();
    //electronQuench::setBirk1(5.1e-3);
    //electronQuench::setEnergyScale(1403.56);
    //electronCerenkov::setA2(-8.75);
    //electronCerenkov::setA3(14.11);
    //electronCerenkov::setA4(0.029);
    //electronResponse::setma(0.99);
    //electronResponse::setmb(7.92e-3);
    //electronResponse::setmc(0);
    //cout << b12->GetChi2() << endl;
    //b12->Plot();

    //for (int i=0; i<8; i++) {
    //    cout << i << " " << electronCerenkov::getCerPE(i) << endl;
    //}

    //electronQuench::setEnergyScale(1800);

    //gammaResponse* gamResArr[9];
    gammaResponse* Cs137 = new gammaResponse("Cs137", 100, 600, 1000);  
    ///gamResArr[0] = Cs137;
    Cs137->LoadData();
    Cs137->preCalculation_onlyNonl();
    std::cout << Cs137->GetName() << " " << Cs137->GetEtrue() << " " << Cs137->GetPEData() << " " << Cs137->GetNonlData() << " " << Cs137->GetChi2_onlyNonl() << endl;

    //////gammaResponse* Mn54 = new gammaResponse("Mn54", 200, 900, 1300);    gamResArr[1] = Mn54;
    //////gammaResponse* Ge68 = new gammaResponse("Ge68", 200, 1100, 1500);   gamResArr[2] = Ge68;
    //////gammaResponse* K40 = new gammaResponse("K40", 200, 1750, 2250);     gamResArr[3] = K40;
    //////gammaResponse* nH = new gammaResponse("nH", 200, 2800, 3500);       gamResArr[4] = nH;
    //////gammaResponse* Co60 = new gammaResponse("Co60", 100, 3200, 3700);   gamResArr[5] = Co60;
    //////gammaResponse* AmBe = new gammaResponse("AmBe", 100, 6000, 6800);   gamResArr[6] = AmBe;
    //////gammaResponse* nC12 = new gammaResponse("nC12", 100, 6700, 7600);   gamResArr[7] = nC12;
    //////gammaResponse* AmC = new gammaResponse("AmC", 100, 8400, 9400);     gamResArr[8] = AmC;
    //////Cs137->SetAmp(250);
    //////Mn54->SetAmp(215);
    //////Ge68->SetAmp(207);
    //////K40->SetAmp(186);
    //////nH->SetAmp(196);
    //////Co60->SetAmp(145);
    //////AmBe->SetAmp(149);
    //////nC12->SetAmp(206);
    //////AmC->SetAmp(158);
    //for (int i=0; i<1; i++) {
    //    gamResArr[i]->LoadData();
    //}

    //for (int i=0; i<1; i++) {
    //    gamResArr[i]->GetChi2();
    //    gamResArr[i]->SaveHist();
    //}


    // water cerenkov for michel electrons
    
    //electronSpectrum* elec500 = new electronSpectrum("500keV", 0.5, 100, 600, 800);
    //elec500->LoadData();
    //elec500->SetAmp(160);
    //elec500->GetChi2();
    //elec500->SaveHist();
    //electronSpectrumFitter *fitter = new electronSpectrumFitter();
    //fitter->LoadData();
    //fitter->GetChiSquare();
    //fitter->Plot(); 

    //michelElectron_water* mic = new michelElectron_water();
    //mic->Initialize();
    //mic->GetChi2();
    //mic->Plot();
    //mic->Compare();

    //michelFitter* mic = new michelFitter();
    //mic->LoadData();
    //mic->GetChiSquare();
    //mic->Plot();

    //michelElectron* mic = new michelElectron();
    //mic->Initialize();
    //cout << mic->GetChi2() << endl;
    //mic->Plot();

    return 1.0;
}
