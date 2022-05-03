#include "gammaResponse.hh"
#include "junoParameters.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"
#include "electronResponse.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <TMath.h>
#include <TFile.h>
#include <TRandom.h>
#include <TF1.h>
#include <TTree.h>
#include <TStopwatch.h>

gammaResponse::gammaResponse(string name, int nBins, double peMin, double peMax) {
    m_name = name;
    m_nBins = nBins;
    m_peMin = peMin;
    m_peMax = peMax;

    //cout << " >>> " << m_name << " " << m_nBins << " bins between [" << m_peMin << ", " << m_peMax <<"]" <<endl;
    hCalc = new TH1D((m_name+"_calc").c_str(), "", m_nBins, m_peMin, m_peMax); 
    hData = new TH1D((m_name+"_data").c_str(), "", m_nBins, m_peMin, m_peMax); 
    func  = new TF1("func", "[0] * TMath::Exp(-(x-[1])*(x-[1])/2/[2]/[2])", m_peMin, m_peMax);

    // init
    m_amp = 250;
    m_loadData  = false;
    m_loadPrm   = false;
    m_doSpecFit = true;


    if (junoParameters::expConfig == "DYB") {
        // DYB Config
        m_Y = 450.87 /2.223;
        m_npeGe68 = 188.6;
        m_sigmaGe68 = 13.92;
    
    }
    
    if (junoParameters::expConfig == "JUNO") {
    // JUNO  Config
    m_npeGe68 = 2*660.8;
    m_sigmaGe68 = 38.41;
    m_Y = 3134.078 /2.223;
    
    }
    
    if (junoParameters::expConfig == "TAO") {
    // JUNO-TAo  Config
    m_npeGe68 = 4.14185e+03;
    m_sigmaGe68 = 7.80287e+01;
    m_Y = 9.89752e+03 /2.223;
    }

    if (junoParameters::expConfig == "Det1") {
        m_npeGe68 = 1.30880e+03;
        m_sigmaGe68 = 3.80672e+01;
        m_Y = 3.11144e+03/2.223;

    }
    if (junoParameters::expConfig == "Det3") {
        m_npeGe68 = 2.61736e+03;
        m_sigmaGe68 = 5.71569e+01;
        m_Y = 6.22116e+03/2.223;

    }
}

gammaResponse::~gammaResponse()
{
    delete hPrmElec;
    delete hPrmPosi;

    delete hCalc;
}


void gammaResponse::LoadData()
{
    cout << " >>> Loading Bared Gamma " << m_name << " Data <<< " << endl;

    ifstream in; in.open(junoParameters::gammaLSNL_File);  // JUNO normal
    //ifstream in; in.open("./data/gamma/gamma_dyb.txt");  // JUNO normal
    if(!in) cout << "Error: No file "<< junoParameters::gammaLSNL_File << std::endl;
    string line;

    double scale = m_Y;
    string tmp_name; double tmp_E, tmp_totPE, tmp_totPESigma, tmp_EvisError, tmp_totPEerr, tmp_totPESigmaerr;
    while(getline(in,line)){
        istringstream ss(line);
        //ss >> tmp_name >> tmp_E >> tmp_totPE >> tmp_totPESigma >> tmp_EvisError ;
        ss >> tmp_name >> tmp_E >> tmp_totPE >> tmp_totPEerr >> tmp_totPESigma >> tmp_totPESigmaerr;
        if(tmp_name == m_name) {
            m_Etrue = tmp_E;
            m_totpeData = tmp_totPE;
            m_nonlData = tmp_totPE/scale/tmp_E;
            //m_nonlDataErr = tmp_totPEerr/scale/tmp_E;
            m_nonlDataErr = 0.001;
            m_Evis = tmp_totPE/scale;
            m_resData = tmp_totPESigma/tmp_totPE;
            m_resDataErr = TMath::Sqrt(tmp_totPESigmaerr*tmp_totPESigmaerr/tmp_totPE/tmp_totPE + tmp_totPEerr*tmp_totPEerr*tmp_totPESigma*tmp_totPESigma/tmp_totPE/tmp_totPE/tmp_totPE/tmp_totPE);
            break;
        }
    }
    in.close();

    if (m_doSpecFit and junoParameters::expConfig == "Det1") {
        TFile* infile  = new TFile(("/junofs/users/miaoyu/simulation/LS_Sim/jobs/Det1/gamma/"+m_name+".root").c_str(), "read") ; 
        int m_totpe;
        TTree* tt = (TTree*)infile->Get("photon");
        tt->SetBranchAddress("totPE", &m_totpe);
        for(int i=0; i<tt->GetEntries(); i++) {
            tt->GetEntry(i);
            hData->Fill(m_totpe);
        }
    }
    else if (m_doSpecFit and junoParameters::expConfig == "Det3") {
        TFile* infile  = new TFile(("/junofs/users/miaoyu/simulation/LS_Sim/jobs/Det3/gamma/"+m_name+".root").c_str(), "read") ; 
        int m_totpe;
        TTree* tt = (TTree*)infile->Get("photon");
        tt->SetBranchAddress("totPE", &m_totpe);
        for(int i=0; i<tt->GetEntries(); i++) {
            tt->GetEntry(i);
            hData->Fill(m_totpe);
        }
    }

    else if (m_doSpecFit and (m_name != "gamma4440" and m_name!="gamma3215" ) ) {
        TFile* infile;
        //TFile* infile = new TFile(("./data/gamma/spectrum/" + m_name + "_totpe.root").c_str(), "read");
        if (junoParameters::expConfig == "JUNO")
            infile = new TFile(("./data/gamma/spectrum/" + m_name + "_totpe_new.root").c_str(), "read");   // JUNO normal
        if (junoParameters::expConfig == "TAO")
            infile = new TFile(("./data/gamma/spectrum/" + m_name + "_totpe_tao.root").c_str(), "read");   // JUNO normal
        if (junoParameters::expConfig == "DYB")
            infile = new TFile(("./data/gamma/spectrum/" + m_name + "_totpe_dyb.root").c_str(), "read");   // JUNO normal
        if(!infile) cout << "No such gamma spectrum file!" << endl;
        int m_totpe;
        //double m_totpe;
        TTree* tt = (TTree*)infile->Get("evt");
        //tt->SetBranchAddress("totalPE", &m_totpe);
        tt->SetBranchAddress("totpe", &m_totpe);
        for(int i=0; i<tt->GetEntries(); i++) {
            tt->GetEntry(i);
            hData->Fill(m_totpe);
        }
    }

    cout << " ****************** Data Size: " << hData->GetEntries() << endl;


    LoadPrmBeta();
}

void gammaResponse::LoadPrmBeta()
{
    //cout << " >>> Load Primary Electron in Single Event for " << m_name << " <<< " << endl;
    string filename = "./data/gamma/" + m_name + "_J19.root";
    TFile* file = new TFile(filename.c_str(), "read");
    if (!file) cout << " No such input file: " << filename << endl;
    hPrmElec = (TH2D*)file->Get((m_name+"_elec").c_str());
    hPrmPosi = (TH2D*)file->Get((m_name+"_posi").c_str());
    
    m_loadPrm = true;
}


void gammaResponse::preCalculation()
{
    //TStopwatch timer;
    //timer.Start();
    //if (not electronResponse::getLoadResol()) electronResponse::loadElecResol();
    //electronResponse::setna(0.94);
    //electronResponse::setnb(0.099);
    //electronResponse::setnc(1.449);
    electronResponse::SetParameters();

    for (int index=0; index<m_nSamples; index++) {
        double tmp_pe = 0;
        double tmp_sigma = 0;
        // electron
        for (int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmElec->GetBinContent(index+1, iSec+1);    
            if (tmp_E == 0) break;
            double single_npe = electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E);
            tmp_pe += single_npe;
            if (junoParameters::pesigmaMode == "kTotal" ) {
                tmp_sigma += TMath::Power(electronResponse::fElecResol->Eval(tmp_E), 2);
                //cout << tmp_E << " " << single_npe << " " << electronResponse::fElecResol->Eval(tmp_E) << " " << electronResponse::fNPESigma->Eval(single_npe) << endl;
            } else if (junoParameters::pesigmaMode == "kNPE") {
                //tmp_sigma += TMath::Power(electronResponse::fNPESigma->Eval(single_npe), 2);             // consider sigma-NPE relationship
                tmp_sigma += TMath::Power(electronResponse::fEvisSigma->Eval(single_npe), 2);             // consider sigma-NPE relationship

            } else if (junoParameters::pesigmaMode == "kSeparate") {

                double sctpe = electronQuench::ScintillatorPE(tmp_E);
                double cerpe = electronCerenkov::getCerPE(tmp_E);
                tmp_sigma += sctpe + electronResponse::fCerPESigma->Eval(cerpe) + 2 * electronResponse::fNtotCov->Eval(sctpe+cerpe);
            }  else if (junoParameters::pesigmaMode == "kNew") {
                tmp_sigma += TMath::Power(electronResponse::fEvisNew->Eval(single_npe), 2);             // consider sigma-NPE relationship
            }

            //tmp_sigma += TMath::Power(electronResponse::gElecResol->Eval(tmp_E), 2);
        }

        // positron
        for(int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmPosi->GetBinContent(index+1, iSec+1);
            if (tmp_E == 0) break;
            double single_npe = electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E);
            tmp_pe += single_npe;
            //tmp_sigma += TMath::Power(electronResponse::fElecResol->Eval(tmp_E), 2);
            if (junoParameters::pesigmaMode == "kTotal" )  {
                tmp_sigma += TMath::Power(electronResponse::fElecResol->Eval(tmp_E), 2);
            } else if (junoParameters::pesigmaMode == "kNPE") {
                //tmp_sigma += TMath::Power(electronResponse::fNPESigma->Eval(single_npe), 2);             // consider sigma-NPE relationship
                tmp_sigma += TMath::Power(electronResponse::fEvisSigma->Eval(single_npe), 2);             // consider sigma-NPE relationship
                //tmp_sigma += TMath::Power(electronResponse::fEvisNew->Eval(single_npe), 2);             // consider sigma-NPE relationship

            } else if (junoParameters::pesigmaMode == "kSeparate") {
                double sctpe = electronQuench::ScintillatorPE(tmp_E);
                double cerpe = electronCerenkov::getCerPE(tmp_E);
                tmp_sigma += sctpe + electronResponse::fCerPESigma->Eval(cerpe) + 2 * electronResponse::fNtotCov->Eval(sctpe+cerpe);
            }  else if (junoParameters::pesigmaMode == "kNew") {
                //tmp_sigma += TMath::Power(electronResponse::fEvisNew1->Eval(single_npe), 2);             // consider sigma-NPE relationship
                tmp_sigma += TMath::Power(electronResponse::fEvisNew->Eval(single_npe), 2);             // consider sigma-NPE relationship
            }
            //tmp_sigma += TMath::Power(electronResponse::gElecResol->Eval(tmp_E), 2);
            tmp_pe += m_npeGe68;
            tmp_sigma += TMath::Power(m_sigmaGe68, 2);

        }

        tmp_sigma = TMath::Sqrt(tmp_sigma);

        m_pemean[index] = tmp_pe;
        m_pesigma[index] = tmp_sigma;

    }

    //timer.Stop();

    //cout << "Total time collapsed during preCalculation : " << timer.RealTime() << " , " << timer.CpuTime() << endl;
}


void gammaResponse::Prediction()
{
    preCalculation();
    double tot_mean = 0;
    for (int i=0; i<m_nSamples; i++)  {
        tot_mean += m_pemean[i];
    }
    tot_mean /= m_nSamples;

    double tot_sigma = 0;
    for (int i=0; i<m_nSamples; i++) {
        tot_sigma += (m_pemean[i]-tot_mean)*(m_pemean[i]-tot_mean) + m_pesigma[i]*m_pesigma[i];
    }
    tot_sigma = TMath::Sqrt(tot_sigma/m_nSamples);

    m_totpeCalc = tot_mean;
    m_totpeSigmaCalc = tot_sigma;

    m_nonlCalc = tot_mean / electronQuench::getEnergyScale() / m_Etrue;
    m_resCalc  = tot_sigma / tot_mean;
}



double gammaResponse::SampleGamEnergy(int index)
{
    if (index >= m_nEvents or index < 0) { cout << "Incorrect Index !" << endl; return 0; } 
    else {
        double tmp_pe = 0;
        double tmp_sigma = 0;
        // electron
        for (int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmElec->GetBinContent(index+1, iSec+1);    
            if (tmp_E == 0) break;
            tmp_pe += electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E);
            tmp_sigma += TMath::Power(electronResponse::gElecResol->Eval(tmp_E), 2);
        }

        // positron
        for(int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmPosi->GetBinContent(index+1, iSec+1);
            if (tmp_E == 0) break;
            tmp_pe += electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E) + 2*660.8;
            tmp_sigma += TMath::Power(electronResponse::gElecResol->Eval(tmp_E), 2);
            tmp_sigma += 2*TMath::Power(27.07, 2);
            
        }

        tmp_sigma = TMath::Sqrt(tmp_sigma);
    
        double sample_pe = gRandom->Gaus(tmp_pe, tmp_sigma);
        return sample_pe;
    }
}




void gammaResponse::calcGamResponse()
{
    if (not electronResponse::getLoadResol()) electronResponse::loadElecResol();
    hCalc->Reset();

    for (int iSample = 0; iSample<m_nSamples; iSample++) {
        int index = gRandom->Integer(5000);
        double sample_pe = SampleGamEnergy(index);
        hCalc->Fill(sample_pe); 
    }
    //hCalc->Fit("gaus", "Q0");
    //double pe_mean  = hCalc->GetFunction("gaus")->GetParameter(1);
    //double pe_sigma = hCalc->GetFunction("gaus")->GetParameter(2);

    //double pe_mean  = hCalc->GetMean();
    //double pe_sigma = hCalc->GetStdDev();

    hCalc->Fit("gaus", "Q0");
    double pe_mean =  hCalc->GetFunction("gaus")->GetParameter(1);
    double pe_sigma = hCalc->GetFunction("gaus")->GetParameter(2);

    m_nonlCalc = pe_mean / electronQuench::getEnergyScale() / m_Etrue;
    m_resCalc = pe_sigma / pe_mean;

}



void gammaResponse::preCalculation_onlyNonl()
{
    //double Y = 3134.078 / 2.223;  // JUNO
    double Y = m_Y;
    double sumpe = 0;
    for (int index=0; index<m_nSamples; index++) {
        double tmp_pe = 0;

        // electron
        for (int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmElec->GetBinContent(index+1, iSec+1);    
            if (tmp_E == 0) break;
            double single_npe = electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E);
            tmp_pe += single_npe;
        }

        // positron
        for(int iSec=0; iSec<100; iSec++) {
            double tmp_E = hPrmPosi->GetBinContent(index+1, iSec+1);
            if (tmp_E == 0) break;
            double single_npe = electronQuench::ScintillatorPE(tmp_E) + electronCerenkov::getCerPE(tmp_E) + m_npeGe68; 
            tmp_pe += single_npe;
        }
        sumpe += tmp_pe;
    }

    m_totpeCalc = sumpe / m_nSamples;
    m_nonlCalc = m_totpeCalc / Y / m_Etrue;
}


double gammaResponse::GetChi2_onlyNonl()
{
    preCalculation_onlyNonl();
    double chi2 = TMath::Power((m_nonlCalc - m_nonlData)/m_nonlDataErr, 2);
    cout << m_name << " " << m_nonlCalc << " " << m_nonlData << " " << m_nonlDataErr << " " << chi2 << endl;
    return chi2;
}




double gammaResponse::GetChi2()
{
    double chi2 = 0;

    //calcGamResponse();
    Prediction();
    
    func->SetParameters(m_amp, m_totpeCalc, m_totpeSigmaCalc);

    if (junoParameters::expConfig == "JUNO" or junoParameters::expConfig=="DYB" or junoParameters::expConfig == "TAO" or junoParameters::expConfig == "Det1" or junoParameters::expConfig=="Det3") {
        if (not m_doSpecFit) 
            chi2 += TMath::Power((m_nonlData - m_nonlCalc)/m_nonlDataErr, 2);

        else{
            for (int i=0; i<m_nBins; i++) {
                double m_err = hData->GetBinError(i);
                if (m_err != 0) {
                    double m_bins = hData->GetBinCenter(i);
                    double m_data = hData->GetBinContent(i);
                    //double m_calc = hCalc->GetBinContent(i);
                    double m_calc = func->Eval(m_bins);
                    //cout << "spectrum fitting " << m_bins << " " << m_data << " " << m_calc << " " << m_err<< endl;

                    chi2 += TMath::Power((m_calc - m_data)/m_err, 2);
                    //cout << hData->GetBinCenter(i) << " " << hCalc->GetBinCenter(i) << " " << m_calc << " " << m_data << " " << m_err << " " << (m_calc- m_data)*(m_calc-m_data)/m_err/m_err<< " "<<chi2 << endl;
                }
            }
        }
    }


    //if (junoParameters::expConfig == "TAO")  {
    //
    //    m_nonlData = hData->GetMean() / m_Etrue / m_Y;
    //    m_nonlDataErr = 0.001;
    //    m_resData = hData->GetStdDev() / hData->GetMean();
    //    m_resDataErr = 0.01;

    //    m_resCalc = m_totpeSigmaCalc / m_totpeCalc;

    //    chi2 += TMath::Power((m_nonlCalc - m_nonlData)/m_nonlDataErr, 2) + TMath::Power((m_resData-m_resCalc)/m_resDataErr, 2);
    //
    //    cout << m_name << " " << m_nonlCalc << " " << m_nonlData << " " << m_resCalc << " " << m_resData << " " << chi2 << endl;
    //}

    //cout << m_name << " " << chi2 << endl;
    return chi2 ;
}


void gammaResponse::SaveHist()
{
    string fileName = m_name + "_"+junoParameters::expConfig+"_hist.root";
    TFile* outfile = new TFile(fileName.c_str(), "recreate");
    hData->Write();
    //hCalc->Write();
    func->Write();
    outfile->Close();
}

