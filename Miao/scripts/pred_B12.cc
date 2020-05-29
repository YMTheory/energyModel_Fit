#include <TMath.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1.h>
#include <TCanvas.h>

TGraph* relative_corr(double *p1, double *p2) {
    TGraph* gRela = new TGraph();
    for(int i=0; i<1000; i++) {
        gRela->SetPoint(i, 13.4/1000*(i+1), p1[i]/p2[i]);
    }
    return gRela;
}

void pred_B12()
{
    double K = 1; // global normalization

    // traditional phase factor 

    // basic parameters:
    double p = 0; // total integral 
    double W = TMath::Sqrt(p*p+1);
    double alpha = 1/137.036;
    double Z = 12; double A = 5;
    double gamma = TMath::Sqrt((1-(alpha*Z)*(alpha*Z)));
    double R = 0.0029*TMath::Power(A,1/3)+0.0063*TMath::Power(A,-1/3)-0.017*TMath::Power(A,-1);
    double E0 = 13.4; //MeV
    double W0 = E0/0.511+1;  //MeV

    // Fermi function part: 
    double F1 = 2*(gamma+1)*TMath::Power(ROOT::Math::tgamma(2*gamma+1),-2)*TMath::Power((2*R),2*(gamma-1));  // energy-independent
    // Table for expansions ...
    double Z0K0L0 = 1; 
    double Z1K1L0 = 3.141593;
    double Z2K0L0 = 0.577216; double Z2K0L1 = -1; double Z2K2L0 = 3.289868;
    double Z3K1L0 = 1.813376; double Z3K1L1 = -3.141593;
    double Z4K0L0 = 0.722126; double Z4K0L1 = -0.827216; double Z4K0L2 = 0.5; double Z4K2L0 = 0.696907; double Z4K2L1 = -3.289868; double Z4K4L0 = -2.164646;
    double Z5K1L0 = 2.268627; double Z5K1L1 = -2.598775; double Z5K1L2 = 1.570786; double Z5K3L0 = -3.776373;
    double Z6K0L0 = 0.730658; double Z6K0L1 = -0.991430; double Z6K0L2 = 0.538608; double Z6K0L3 = -0.166667; 
    double Z6K2L0 = 0.569598; double Z6K2L1 = -1.519374; double Z6K2L2 = 1.644934; double Z6K4L0 = -4.167149; double Z6K4L1 = 2.164646; double Z6K6L0 = 2.034686;
    double Z7K1L0 = 2.295429; double Z7K1L1 = -3.114670; double Z7K1L2 = 1.692086; double Z7K1L3 = -0.523599; double Z7K3L0 = -5.674039; double Z7K3L1 = 3.776373; double Z7K5L0=3.257650;
    double Z8K0L0 = 0.752192; double Z8K0L1 = -1.061466; double Z8K0L2 = 0.661617; double Z8K0L3 = -0.221203; double Z8K0L4 = 0.041667;
    double Z8K2L0 = -0.180873; double Z8K2L1 = -1.155058; double Z8K2L2 = 1.170920; double Z8K2L3 = -0.548311;
    double Z8K4L0 = -4.653076; double Z8K4L1 = 4.708310; double Z8K4L2 = -1.082323; double Z8K6L0=6.179487; double Z8K6L1 = -2.034686; double Z8K8L0 = -2.008155;
    double Z9K1L0 = 2.363082; double Z9K1L1 = -3.334694; double Z9K1L2 = 2.078531; double Z9K1L3 = -0.694928; double Z9K1L4 = 0.130900;
    double Z9K3L0 = -8.119890; double Z9K3L1 = 6.618132; double Z9K3L2 = -1.888187; 
    double Z9K5L0 = 8.959546; double Z9K5L1 = -3.257605; double Z9K7L0 = -3.167823;
    double Z10K0L0 = 0.768082; double Z10K0L1 = -1.124905; double Z10K0L2 = -0.745425; double Z10K0L3 = -0.286256; double Z10K0L4 = 0.065717; double Z10K0L5 = -0.008333;
    double Z10K2L0 = -0.881488; double Z10K2L1 = -0.305660; double Z10K2L2 = 0.973067; double Z10K2L3 = -0.527385; double Z10K2L4 = 0.137078;
    double Z10K4L0 = -4.672377; double Z10K4L1 = 5.965444; double Z10K4L2 = -2.624736; double Z10K4L3 = 0.360774;
    double Z10K6L0 = 10.923585; double Z10K6L1 = -6.688159; double Z10K6L2 = 1.017343; 
    double Z10K8L0 = -8.164857; double Z10K8L1 = 2.008155; double Z10K10L0 = 2.001989;
    double Z11K1L0 = 2.413001; double Z11K1L1 = -3.533993; double Z11K1L2 = 2.341823; double Z11K1L3 = -0.899301; double Z11K1L4 = 0.206457; double Z11K1L5=-0.026180;
    double Z11K3L0 = -10.543505; double Z11K3L1 = 10.010446; double Z11K3L2 = -3.781113; double Z11K3L3 = 0.629396;
    double Z11K5L0 = 17.003454; double Z11K5L1 = -9.773947;double Z11K5L2= 1.628802; double Z11K7L0=-12.056502; double Z11K7L1 = 3.167823; double Z11K9L0=3.147902;

    double KE[1000]; double spec[1000]; double fm_spec[1000]; double fsc_spec[1000]; double wi_spec[1000]; double se_spec[1000]; double wm_spec[1000];
    double scale1 = 0;

    for(int i=0; i<1000; i++) {
        double betaE = E0/1000.*(i+1);
        double betaW = betaE/0.511+1;
        double betaP = TMath::Sqrt(betaW*betaW-1);
        KE[i] = betaE;

        // traditional phase space : 
        double sp = betaP*betaP*(W0-betaW)*(W0-betaW);

        double t1 = alpha*Z; double t2 = betaW/betaP; double t3=TMath::Log(betaP);

        // long long long expanesion ...
        double F2 = Z0K0L0*TMath::Power(t1,0)*TMath::Power(t2,0)*TMath::Power(t3,0)
                + Z1K1L0*TMath::Power(t1,1)*TMath::Power(t2,1)*TMath::Power(t3,0)
                + Z2K0L0*TMath::Power(t1,2)*TMath::Power(t2,0)*TMath::Power(t3,0)
                + Z2K0L1*TMath::Power(t1,2)*TMath::Power(t2,0)*TMath::Power(t3,1)
                + Z2K2L0*TMath::Power(t1,2)*TMath::Power(t2,2)*TMath::Power(t3,0)
                + Z3K1L0*TMath::Power(t1,3)*TMath::Power(t2,1)*TMath::Power(t3,0)
                + Z3K1L1*TMath::Power(t1,3)*TMath::Power(t2,1)*TMath::Power(t3,1)
                + Z4K0L0*TMath::Power(t1,4)*TMath::Power(t2,0)*TMath::Power(t3,0)
                + Z4K0L1*TMath::Power(t1,4)*TMath::Power(t2,0)*TMath::Power(t3,1)
                + Z4K0L2*TMath::Power(t1,4)*TMath::Power(t2,0)*TMath::Power(t3,2)
                + Z4K2L0*TMath::Power(t1,4)*TMath::Power(t2,2)*TMath::Power(t3,0)
                + Z4K2L1*TMath::Power(t1,4)*TMath::Power(t2,2)*TMath::Power(t3,1)
                + Z4K4L0*TMath::Power(t1,4)*TMath::Power(t2,4)*TMath::Power(t3,0)
                + Z5K1L0*TMath::Power(t1,5)*TMath::Power(t2,1)*TMath::Power(t3,0)
                + Z5K1L1*TMath::Power(t1,5)*TMath::Power(t2,1)*TMath::Power(t3,1)
                + Z5K1L2*TMath::Power(t1,5)*TMath::Power(t2,1)*TMath::Power(t3,2)
                + Z5K3L0*TMath::Power(t1,5)*TMath::Power(t2,3)*TMath::Power(t3,0)
                + Z6K0L0*TMath::Power(t1,6)*TMath::Power(t2,0)*TMath::Power(t3,0)
                + Z6K0L1*TMath::Power(t1,6)*TMath::Power(t2,0)*TMath::Power(t3,1)
                + Z6K0L2*TMath::Power(t1,6)*TMath::Power(t2,0)*TMath::Power(t3,2)
                + Z6K0L3*TMath::Power(t1,6)*TMath::Power(t2,0)*TMath::Power(t3,3)
                + Z6K2L0*TMath::Power(t1,6)*TMath::Power(t2,2)*TMath::Power(t3,0)
                + Z6K2L1*TMath::Power(t1,6)*TMath::Power(t2,2)*TMath::Power(t3,1)
                + Z6K2L2*TMath::Power(t1,6)*TMath::Power(t2,2)*TMath::Power(t3,2)
                + Z6K4L0*TMath::Power(t1,6)*TMath::Power(t2,4)*TMath::Power(t3,0)
                + Z6K4L1*TMath::Power(t1,6)*TMath::Power(t2,4)*TMath::Power(t3,1)
                + Z6K6L0*TMath::Power(t1,6)*TMath::Power(t2,6)*TMath::Power(t3,0);
        //        + Z7K1L0*TMath::Power(t1,7)*TMath::Power(t2,1)*TMath::Power(t3,0)
        //        + Z7K1L1*TMath::Power(t1,7)*TMath::Power(t2,1)*TMath::Power(t3,1)
        //        + Z7K1L2*TMath::Power(t1,7)*TMath::Power(t2,1)*TMath::Power(t3,2)
        //        + Z7K1L3*TMath::Power(t1,7)*TMath::Power(t2,1)*TMath::Power(t3,3)
        //        + Z7K3L0*TMath::Power(t1,7)*TMath::Power(t2,3)*TMath::Power(t3,0)
        //        + Z7K3L1*TMath::Power(t1,7)*TMath::Power(t2,3)*TMath::Power(t3,1)
        //        + Z7K5L0*TMath::Power(t1,7)*TMath::Power(t2,5)*TMath::Power(t3,0)
        //        + Z8K0L0*TMath::Power(t1,8)*TMath::Power(t2,0)*TMath::Power(t3,0)
        //        + Z8K0L1*TMath::Power(t1,8)*TMath::Power(t2,0)*TMath::Power(t3,1)
        //        + Z8K0L2*TMath::Power(t1,8)*TMath::Power(t2,0)*TMath::Power(t3,2)
        //        + Z8K0L3*TMath::Power(t1,8)*TMath::Power(t2,0)*TMath::Power(t3,3)
        //        + Z8K0L4*TMath::Power(t1,8)*TMath::Power(t2,0)*TMath::Power(t3,4)
        //        + Z8K2L0*TMath::Power(t1,8)*TMath::Power(t2,2)*TMath::Power(t3,0)
        //        + Z8K2L1*TMath::Power(t1,8)*TMath::Power(t2,2)*TMath::Power(t3,1)
        //        + Z8K2L2*TMath::Power(t1,8)*TMath::Power(t2,2)*TMath::Power(t3,2)
        //        + Z8K2L3*TMath::Power(t1,8)*TMath::Power(t2,2)*TMath::Power(t3,3)
        //        + Z8K4L0*TMath::Power(t1,8)*TMath::Power(t2,4)*TMath::Power(t3,0)
        //        + Z8K4L1*TMath::Power(t1,8)*TMath::Power(t2,4)*TMath::Power(t3,1)
        //        + Z8K4L2*TMath::Power(t1,8)*TMath::Power(t2,4)*TMath::Power(t3,2)
        //        + Z8K6L0*TMath::Power(t1,8)*TMath::Power(t2,6)*TMath::Power(t3,0)
        //        + Z8K6L1*TMath::Power(t1,8)*TMath::Power(t2,6)*TMath::Power(t3,1)
        //        + Z8K8L0*TMath::Power(t1,8)*TMath::Power(t2,8)*TMath::Power(t3,0)
        //        + Z9K1L0*TMath::Power(t1,9)*TMath::Power(t2,1)*TMath::Power(t3,0)
        //        + Z9K1L1*TMath::Power(t1,9)*TMath::Power(t2,1)*TMath::Power(t3,1)
        //        + Z9K1L2*TMath::Power(t1,9)*TMath::Power(t2,1)*TMath::Power(t3,2)
        //        + Z9K1L3*TMath::Power(t1,9)*TMath::Power(t2,1)*TMath::Power(t3,3)
        //        + Z9K1L4*TMath::Power(t1,9)*TMath::Power(t2,1)*TMath::Power(t3,4)
        //        + Z9K3L0*TMath::Power(t1,9)*TMath::Power(t2,3)*TMath::Power(t3,0)
        //        + Z9K3L1*TMath::Power(t1,9)*TMath::Power(t2,3)*TMath::Power(t3,1)
        //        + Z9K3L2*TMath::Power(t1,9)*TMath::Power(t2,3)*TMath::Power(t3,2)
        //        + Z9K5L0*TMath::Power(t1,9)*TMath::Power(t2,5)*TMath::Power(t3,0)
        //        + Z9K5L1*TMath::Power(t1,9)*TMath::Power(t2,5)*TMath::Power(t3,1)
        //        + Z9K7L0*TMath::Power(t1,9)*TMath::Power(t2,7)*TMath::Power(t3,0)
        //        + Z10K0L0*TMath::Power(t1,10)*TMath::Power(t2,0)*TMath::Power(t3,0)
        //        + Z10K0L1*TMath::Power(t1,10)*TMath::Power(t2,0)*TMath::Power(t3,1)
        //        + Z10K0L2*TMath::Power(t1,10)*TMath::Power(t2,0)*TMath::Power(t3,2)
        //        + Z10K0L3*TMath::Power(t1,10)*TMath::Power(t2,0)*TMath::Power(t3,3)
        //        + Z10K0L4*TMath::Power(t1,10)*TMath::Power(t2,0)*TMath::Power(t3,4)
        //        + Z10K0L5*TMath::Power(t1,10)*TMath::Power(t2,0)*TMath::Power(t3,5)
        //        + Z10K2L0*TMath::Power(t1,10)*TMath::Power(t2,2)*TMath::Power(t3,0)
        //        + Z10K2L1*TMath::Power(t1,10)*TMath::Power(t2,2)*TMath::Power(t3,1)
        //        + Z10K2L2*TMath::Power(t1,10)*TMath::Power(t2,2)*TMath::Power(t3,2)
        //        + Z10K2L3*TMath::Power(t1,10)*TMath::Power(t2,2)*TMath::Power(t3,3)
        //        + Z10K2L4*TMath::Power(t1,10)*TMath::Power(t2,2)*TMath::Power(t3,4)
        //        + Z10K4L0*TMath::Power(t1,10)*TMath::Power(t2,4)*TMath::Power(t3,0)
        //        + Z10K4L1*TMath::Power(t1,10)*TMath::Power(t2,4)*TMath::Power(t3,1)
        //        + Z10K4L2*TMath::Power(t1,10)*TMath::Power(t2,4)*TMath::Power(t3,2)
        //        + Z10K4L3*TMath::Power(t1,10)*TMath::Power(t2,4)*TMath::Power(t3,3)
        //        + Z10K6L0*TMath::Power(t1,10)*TMath::Power(t2,6)*TMath::Power(t3,0)
        //        + Z10K6L1*TMath::Power(t1,10)*TMath::Power(t2,6)*TMath::Power(t3,1)
        //        + Z10K6L2*TMath::Power(t1,10)*TMath::Power(t2,6)*TMath::Power(t3,2)
        //        + Z10K8L0*TMath::Power(t1,10)*TMath::Power(t2,8)*TMath::Power(t3,0)
        //        + Z10K8L1*TMath::Power(t1,10)*TMath::Power(t2,8)*TMath::Power(t3,1);
        //        //+ Z10K10L0*TMath::Power(t1,10)*TMath::Power(t2,10)*TMath::Power(t3,0);
        //        //+ Z11K1L0*TMath::Power(t1,11)*TMath::Power(t2,1)*TMath::Power(t3,0)
        //        //+ Z11K1L1*TMath::Power(t1,11)*TMath::Power(t2,1)*TMath::Power(t3,1)
        //        //+ Z11K1L2*TMath::Power(t1,11)*TMath::Power(t2,1)*TMath::Power(t3,2)
        //        //+ Z11K1L3*TMath::Power(t1,11)*TMath::Power(t2,1)*TMath::Power(t3,3)
        //        //+ Z11K1L4*TMath::Power(t1,11)*TMath::Power(t2,1)*TMath::Power(t3,4)
        //        //+ Z11K1L5*TMath::Power(t1,11)*TMath::Power(t2,1)*TMath::Power(t3,5)
        //        //+ Z11K3L0*TMath::Power(t1,11)*TMath::Power(t2,3)*TMath::Power(t3,0)
        //        //+ Z11K3L1*TMath::Power(t1,11)*TMath::Power(t2,3)*TMath::Power(t3,1)
        //        //+ Z11K3L2*TMath::Power(t1,11)*TMath::Power(t2,3)*TMath::Power(t3,2)
        //        //+ Z11K3L3*TMath::Power(t1,11)*TMath::Power(t2,3)*TMath::Power(t3,3)
        //        //+ Z11K5L0*TMath::Power(t1,11)*TMath::Power(t2,5)*TMath::Power(t3,0)
        //        //+ Z11K5L1*TMath::Power(t1,11)*TMath::Power(t2,5)*TMath::Power(t3,1)
        //        //+ Z11K5L2*TMath::Power(t1,11)*TMath::Power(t2,5)*TMath::Power(t3,2)
        //        //+ Z11K7L0*TMath::Power(t1,11)*TMath::Power(t2,7)*TMath::Power(t3,0)
        //        //+ Z11K7L1*TMath::Power(t1,11)*TMath::Power(t2,7)*TMath::Power(t3,1)
        //        //+ Z11K9L0*TMath::Power(t1,11)*TMath::Power(t2,9)*TMath::Power(t3,0);
        //        
    
        // Finiti-size corrections : 
        double L0 = 1+13*(alpha*Z)*(alpha*Z)/60 - betaW*R*alpha*Z*(41-26*gamma)/(15*(2*gamma-1)) - alpha*Z*R*gamma*(17-2*gamma)/(30*betaW*(2*gamma-1));

        // weak-interaction finite-size effects:
        double C0 = -233/630.*(alpha*Z)*(alpha*Z) - (W0*R)*(W0*R)/5 + 2./35*W0*R*alpha*Z;
        double C1 = -21./35*R*alpha*Z + 4./9*W0*R*R;
        double C2 = -4./9*R*R;
        double C = 1 + C0 +C1*betaW + C2*betaW*betaW;

        // screening correciton 
        double Z_bar = Z - 1;
        double V0 = alpha*alpha*TMath::Power(Z_bar, 0.75)*1.4584000;  // linear interpolation from table
        double betaW_bar = betaW-V0;
        double betaP_bar = TMath::Sqrt(betaW_bar*betaW_bar-1);
        double y = alpha*Z*betaW/betaP;
        double y_bar = alpha*Z*betaW_bar/betaP_bar;
        double S1 = betaW_bar/betaW * TMath::Power(betaP_bar/betaP, 2*gamma-1)*TMath::Exp(TMath::Pi()*(y_bar-y))*TMath::Power(ROOT::Math::tgamma(2*gamma+1),-2);
        double tt1 = alpha*Z; double tt2 = betaW_bar/betaP_bar;
        // expansion table ...
        double Z0K0 = 1;
        double Z2K0 = 0.577216; double Z2K2 = -1.644934;
        double Z4K0 = 0.722126; double Z4K2 = -2.151539; double Z4K4 = 1.894066;
        double Z6K0 = 0.730658; double Z6K2 = -2.993953; double Z6K4 = 4.107516; double Z6K6 = -1.971102;
        double S2 = Z0K0*TMath::Power(tt1, 0)*TMath::Power(tt2, 0)
            + Z2K0*TMath::Power(tt1,2)*TMath::Power(tt2,0)
            + Z2K2*TMath::Power(tt1,2)*TMath::Power(tt2,2)
            + Z4K0*TMath::Power(tt1,4)*TMath::Power(tt2,0)
            + Z4K2*TMath::Power(tt1,4)*TMath::Power(tt2,2)
            + Z4K4*TMath::Power(tt1,4)*TMath::Power(tt2,4)
            + Z6K0*TMath::Power(tt1,6)*TMath::Power(tt2,0)
            + Z6K2*TMath::Power(tt1,6)*TMath::Power(tt2,2)
            + Z6K4*TMath::Power(tt1,6)*TMath::Power(tt2,4)
            + Z6K6*TMath::Power(tt1,6)*TMath::Power(tt2,6);
        double S = 1;
        if(betaW>V0) { S = S1*S2; } 
 
        // weak magnetism :  1+delta_WM * W
        // it is a linear correction with a coefficient of energy spectra
        double deltaWM = 0.0048;

        spec[i] = sp;
        fm_spec[i] = sp*F1*F2;
        fsc_spec[i] = sp*F1*F2*L0;
        wi_spec[i] = sp*F1*F2*L0*C;
        se_spec[i] = sp*F1*F2*L0*C*S;
        wm_spec[i] = se_spec[i]*(1+deltaWM);

        scale1 += wm_spec[i];
        
    }

    // scale all spectra 
    double s1 = spec[500];
    double s2 = fm_spec[500];
    double s3 = fsc_spec[500];
    double s4 = wi_spec[500];
    double s5 = se_spec[500];
    double s6 = wm_spec[500];
    for(int i=0; i<1000; i++) {
        fm_spec[i] = fm_spec[i]/s2*s1;
        fsc_spec[i] = fsc_spec[i]/s3*s1;
        wi_spec[i] = wi_spec[i]/s4*s1;
        se_spec[i] = se_spec[i]/s5*s1;
        wm_spec[i] = wm_spec[i]/s6*s1;
    }

    TCanvas* c1 = new TCanvas(); c1->cd();
    TGraph* gFmCorr = relative_corr(fm_spec, spec);
    TGraph* gFsCorr = relative_corr(fsc_spec, spec);
    TGraph* gWiCorr = relative_corr(wi_spec, spec);
    TGraph* gSeCorr = relative_corr(se_spec, spec);
    TGraph* gWmCorr = relative_corr(wm_spec, spec);

    gFmCorr->SetLineColor(kCyan+1);
    gFsCorr->SetLineColor(kRed);
    gWiCorr->SetLineColor(kOrange+1);
    gSeCorr->SetLineColor(kGreen+1);
    gWmCorr->SetLineColor(kPink+1);

    gFmCorr->SetTitle("B12 Spectrum; beta kinetic energy/MeV; a.u");
    gFmCorr->Draw("AL");
    gFsCorr->Draw("L SAME");
    gWiCorr->Draw("L SAME");
    gSeCorr->Draw("L SAME");
    gWmCorr->Draw("L SAME");

    TLegend* ll = new TLegend();
    ll->AddEntry(gFmCorr, "w/ Fermi Function", "l");
    ll->AddEntry(gFsCorr, "w/ Finite-size Corrections", "l");
    ll->AddEntry(gWiCorr, "w/ weak-interaction finite size corrections", "l");
    ll->AddEntry(gSeCorr, "w/ screening corrections", "l");
    ll->AddEntry(gWmCorr, "w/ weak magnetism corrections", "l");
    ll->Draw("SAME");


    TCanvas* c2 = new TCanvas(); c2->cd();

    TGraph* gSpec = new TGraph(1000, KE, spec);
    TGraph* gFmSpec = new TGraph(1000, KE, fm_spec);
    TGraph* gFscSpec = new TGraph(1000, KE, fsc_spec);
    TGraph* gWiSpec = new TGraph(1000, KE, wi_spec);
    TGraph* gSeSpec = new TGraph(1000, KE, se_spec);
    TGraph* gWmSpec = new TGraph(1000, KE, wm_spec);

    gFmSpec->SetLineColor(kCyan+1);
    gSpec->SetLineColor(kBlue);
    gFscSpec->SetLineColor(kRed);
    gWiSpec->SetLineColor(kOrange+1);
    gSeSpec->SetLineColor(kGreen+1);
    gWmSpec->SetLineColor(kPink+1);

    gFmSpec->SetTitle("B12 Spectrum; beta kinetic energy/MeV; a.u");
    gFmSpec->Draw("AL");
    gFscSpec->Draw("L SAME");
    gSpec->Draw("L SAME");
    gWiSpec->Draw("L SAME");
    gSeSpec->Draw("L SAME");
    gWmSpec->Draw("L SAME");

    TLegend* led = new TLegend();
    led->AddEntry(gSpec, "phase space facto", "l");
    led->AddEntry(gFmSpec, "w/ Fermi Function", "l");
    led->AddEntry(gFscSpec, "w/ Finite-size Corrections", "l");
    led->AddEntry(gWiSpec, "w/ weak-interaction finite size corrections", "l");
    led->AddEntry(gSeSpec, "w/ screening corrections", "l");
    led->AddEntry(gWmSpec, "w/ weak magnetism corrections", "l");
    led->Draw("SAME");

    
    TH1D* hSimul = new TH1D("Simul", "", 200, 0, E0);
    ifstream in;
    in.open("./B12_edep.txt");
    string line; double tmp;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmp;
        hSimul->Fill(tmp);
    }
    double scale2 = hSimul->GetEntries();
    hSimul->Scale(scale1/scale2/5);
    hSimul->SetLineColor(kMagenta);
    hSimul->Draw("SAME");

}
