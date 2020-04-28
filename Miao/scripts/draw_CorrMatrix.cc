void draw_CorrMatrix()
{
    gStyle->SetOptStat(0);

    TH2F* CorrMatrix = new TH2F("CorrMatrix", "", 3, 0, 3, 3, 0, 3);
    CorrMatrix->SetBinContent(1,3,1.000);
    CorrMatrix->SetBinContent(2,2,1.000);
    CorrMatrix->SetBinContent(3,1,1.000);
    CorrMatrix->SetBinContent(1,2,0.953);
    CorrMatrix->SetBinContent(2,3,0.953);
    CorrMatrix->SetBinContent(1,1,-0.933);
    CorrMatrix->SetBinContent(3,3,-0.933);
    CorrMatrix->SetBinContent(2,1,-0.824);
    CorrMatrix->SetBinContent(3,2,-0.824);
    CorrMatrix->Draw("TEXT COLZ");

}
