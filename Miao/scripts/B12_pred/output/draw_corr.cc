#include <TGraph.h>
#include <TFile.h>
#include <TLegend.h>

using namespace std;

void draw_corr()
{
    TGraph* gWeakFiniteSize;
    TGraph* gElecFiniteSize;
    TGraph* gScreening;
    TGraph* gWeakMag;

    TFile* file1 = TFile::Open("./FiniteSize.root");
    gElecFiniteSize = (TGraph*)file1->Get("ElecFiniteSize");
    gWeakFiniteSize = (TGraph*)file1->Get("WeakFiniteSize");

    TFile* file2 = TFile::Open("./screening.root");
    gScreening = (TGraph*)file2->Get("Screening");

    TFile* file3 = TFile::Open("./WeakMagnetism.root");
    gWeakMag = (TGraph*)file3->Get("WeakMagnetism");

    gElecFiniteSize->SetLineColor(kGreen+2);
    gWeakFiniteSize->SetLineColor(kOrange+1);
    gScreening->SetLineColor(kRed+1);
    gWeakMag->SetLineColor(kBlue+1);
    gElecFiniteSize->SetLineWidth(2);
    gWeakFiniteSize->SetLineWidth(2);
    gScreening->SetLineWidth(2);
    gWeakMag->SetLineWidth(2);

    gElecFiniteSize->GetYaxis()->SetRangeUser(0.9,1.4);

    gElecFiniteSize->Draw("AL");
    gWeakFiniteSize->Draw("L SAME");
    gScreening->Draw("L SAME");
    gWeakMag->Draw("L SAME");

    TLegend* led = new TLegend();
    led->AddEntry(gElecFiniteSize, "electromagnetic finite size", "l");
    led->AddEntry(gWeakFiniteSize, "weak finite size", "l");
    led->AddEntry(gWeakMag, "weak magnetism", "l");
    led->AddEntry(gScreening, "screening effect", "l");
    led->Draw("SAME");
}