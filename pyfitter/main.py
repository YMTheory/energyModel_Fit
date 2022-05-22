from electronResponse import electronResponse
from gammaResponse import gamma
from fitter import fitter

if __name__ == "__main__" :
    #electronResponse.initialize()
    #E, kB, Ysct, kC = 1, 6.5e-3, 1400, 1.0
    #print(E, electronResponse.get_Nsct(E, kB, Ysct), electronResponse.get_Ncer(E, kC), electronResponse.get_Nsigma(E))

    #Cs137 = gamma("Cs137", 0.662)
    #Cs137.load_npe_spectrum()
    #Cs137.load_gamma_samples()
    #Cs137.prediction()

    myFitter = fitter()
    myFitter.fit()
    myFitter.plot()
