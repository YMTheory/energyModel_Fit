#include <iostream>
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"

#include "junoEnergyModel.hh"
#include "junoParameters.hh"
#include "junoGlobalFit.hh"
#include "junoGammaData.hh"
#include "junoGammaPeak.hh"

using namespace std;

void SetStyle()
{
	TStyle *style = new TStyle("Modern","Modern Style");
	style->SetTitleFont(43,"xyz");
	style->SetLabelFont(43,"xyz");
	style->SetLegendFont(43);
	style->SetLabelSize(19,"xyz");
	style->SetTitleSize(21,"xyz");
	style->SetLegendBorderSize(0);
	style->SetLegendFillColor(kRed);
	style->SetStatStyle(0);
	style->SetLineColor(kBlue+1);
	style->SetMarkerStyle(20);
	gROOT->SetStyle("Modern"); 
}
