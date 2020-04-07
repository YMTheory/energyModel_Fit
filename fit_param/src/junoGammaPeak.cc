#include "junoGammaPeak.hh"

int    junoGammaPeak::s_count    = 0;

junoGammaPeak::junoGammaPeak ()  {  s_count++; }

junoGammaPeak::junoGammaPeak( string peakName,
                              string pdfName,
                              double eTru_total,
                              double eTru_single)
{
    s_count++;
    Init(peakName,pdfName,eTru_total,eTru_single);
}

junoGammaPeak::~junoGammaPeak()  {  s_count--; }

void junoGammaPeak::Init( string peakName,
                        string pdfName,
                        double eTru_total,
                        double eTru_single)
{
    m_name           =  peakName;
    m_eTru_total     =  eTru_total;
    m_eTru_single    =  eTru_single;
    m_eVisError      = 0;

    TFile file(junoParameters::gammaPdf_file.c_str(),"read");
    TGraph* gGammaPdf = (TGraph*)file.Get(pdfName.c_str());
    double* tmp_eTru = gGammaPdf->GetX();
    double* tmp_prob = gGammaPdf->GetY();
    for(int i=0; i<gGammaPdf->GetN(); i++)  {
        m_pdf_eTru[i] = tmp_eTru[i];
        m_pdf_prob[i] = tmp_prob[i];
        if(m_pdf_prob[i] > 0)  m_nPdf = i;
    }

}

void junoGammaPeak::SetEVis (double val) { 
    m_eVis = val;
    m_dataScintNL = m_eVis / m_eTru_total;
}


void junoGammaPeak::SetEVisError (double val) {
    m_eVisError = val;
}

void junoGammaPeak::UpdateTheoNL () {
    
    m_eVis = 0;
    double eTru1 = 0;
    double prob1 = 0;
    double eTru2 = 0;
    double prob2 = 0;
    double sum  = 0;
    for(int i=0; i<m_nPdf-1; i++) {
        eTru1    = m_pdf_eTru[i];
        prob1    = m_pdf_prob[i];
        eTru2    = m_pdf_eTru[i+1];
        prob2    = m_pdf_prob[i+1];
        sum    += (eTru1*prob1 + eTru2*prob2) * (eTru2-eTru1)/2.;
        m_eVis += (eTru1*prob1*junoEnergyModel::ScintillatorNL(eTru1) + eTru2*prob2*junoEnergyModel::ScintillatorNL(eTru2))*(eTru2-eTru1)/2.;
    }
    m_eVisError = 0.002;
    //cout << "Final: evis " << m_eVis << "  sum " << sum << endl;
    m_eVis *= m_eTru_total/sum;
    //std::cout << "eVis : " << m_eVis << "  eTru_total: " << m_eTru_total << endl;
    m_theoScintNL = m_eVis/m_eTru_total;
}


void junoGammaPeak::UpdateDataNL ()  {
    UpdateTheoNL ();
    for( int i=0; i<8000; i++ )  {
        double eTru = i * 0.001;
        double eVis = eTru * junoEnergyModel::ScintillatorNL(eTru);
        if( eVis>=m_eVis )
        {
            //std::cout << " found " << eVis << std::endl;
            m_dataScintNL = eVis/m_eTru_total;
            return;
        }
    }
}

double junoGammaPeak::GetChi2()
{
    UpdateTheoNL();
    double error2 = m_eVisError*m_eVisError;
    double delta = m_theoScintNL - m_dataScintNL;
    //std::cout << " =============== " << std::endl;
    //std::cout << "theoScintNL: " << m_theoScintNL << "  dataScintNL: " << m_dataScintNL <<std::endl;
    //std::cout << "Chi2 : " << pow(delta, 2)/error2;
    //std::cout << " =============== " << std::endl;
    return pow( delta, 2 ) / error2;

}

