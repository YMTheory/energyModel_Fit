void response_function()
{
    double npe   = 3134.078;
    double sigma = 71.101;

    double alpha = (sigma*sigma+npe*npe) / (TMath::Power(sigma, 4)*(2+3/npe) + 4*npe*npe*(sigma*sigma-2) + 2*npe*(6*sigma*sigma-1));
    double beta  = alpha / npe/npe;

    cout << alpha << " " << beta << endl;

    double gamma = 2*TMath::Power(beta, alpha) / TMath::tgamma(alpha) *  
}
