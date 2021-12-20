import numpy as np
import matplotlib.pyplot as plt
import uproot as up
import ROOT

def theory(E, A):
    Emax = 5.21563e+01
    rho = 7.30525e-01
    x = E/Emax

    return A*0.5*(12*x*x-12*x*x*x + rho * (32/3.*x*x*x -8*x*x) );

def theory1(E):
    A = 100
    mu = 105.6583745
    me = 0.511
    x = 2*E/mu
    alpha = 1/137.035999173
    powterm = (2*alpha/np.pi) * (np.log(mu*x/me)-1)
    return A * x**2 * (6-4*x) * pow((1-x)/(3*x/7), powterm)



def genTheoHist():
    hh = ROOT.TH1F("hh0", "", 7000, 0, 70)
    for i in range(7000):
        Etrue = 0.005 + 0.01*i
        if Etrue > 105.6583745 / 2:
            hh.SetBinContent(i+1, 0)
            continue
        hh.SetBinContent(i+1, theory1(Etrue))

    ff = ROOT.TFile("michel_theo.root", "recreate")
    hh.Write()
    ff.Close()




ratio = 1.53921e+01
p1 = 1.02133e+00 
p2 = 1.26506e-03
#ratio = 15.7
#p1 = 1.76105 
#p2 = 0.00197229

m_peMin = 0
m_peMax = 1540


def gaus(x, mu, sigma):
    return 1/(np.sqrt(2*np.pi)*sigma) * np.exp( -(x-mu)**2/2/sigma**2)



def loadTheo():
    iff = ROOT.TFile("/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/data/electron/michel/Michel_hist.root", "read")
    hh0 = iff.Get("hh0")
    
    if(0):
        bince, binco = [], []
        for i in range(hh0.GetNbinsX()):
            bince.append(hh0.GetBinCenter(i+1))
            binco.append(hh0.GetBinContent(i+1))
        plt.plot(bince, binco, "-")
        plt.show()

    print(hh0)
    return hh0



def loadData():
    totpe = up.open("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/water/michel/michel_totpe_water1.root")["michel"]["totpe"].array()
    hist = ROOT.TH1D("data", "", 100, m_peMin, m_peMax)
    for j in totpe:
        hist.Fill(j)
    hist.Scale(1./hist.GetEntries())
    bince, binco = [], []
    for i in range(100):
        bince.append(hist.GetBinCenter(i+1))
        binco.append(hist.GetBinContent(i+1))

    return bince, binco



def applyCerNPE():
    
    iff = ROOT.TFile("/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/data/electron/michel/Michel_hist.root", "read")
    hist = iff.Get("hh0")

    m_nBin = hist.GetNbinsX()
    m_binWidth = hist.GetBinCenter(2) - hist.GetBinCenter(1)
    print("fine binning width = %.3f MeV/bin with total %d bins ." %(m_binWidth, m_nBin) )

    m_peBinWidth = ( m_peMax - m_peMin ) / m_nBin

    m_eVis = [0 for i in range(m_nBin)]

    for i in range(m_nBin):
        eTru = hist.GetBinCenter(i+1)
        value = hist.GetBinContent(i+1)
        tmp_pe = eTru * ratio
        tmp_sigma = np.sqrt(p1*tmp_pe + p2*tmp_pe**2)
        minBin = int((tmp_pe - 5*tmp_sigma)/m_peBinWidth)
        maxBin = int((tmp_pe + 5*tmp_sigma)/m_peBinWidth)
        if minBin < 0:
            minBin = 0
        if maxBin > m_nBin-1:
            maxBin = m_nBin-1
        #print(eTru, value, tmp_pe, tmp_sigma, tmp_pe-5*tmp_sigma, minBin, tmp_pe+5*tmp_sigma, maxBin)


        for j in range(minBin, maxBin+1, 1):
            tmp_center = m_binWidth/2 + m_peBinWidth * j
            prob = gaus(tmp_center, tmp_pe, tmp_sigma)

            m_eVis[j] += prob * value / 6


    binCenter = []
    for i in range(m_nBin):
        binCenter.append( i*m_peBinWidth + m_peBinWidth/2.)

    databince, databinco = loadData()
    #m_eTru = []
    #for i in range(m_nBin):
    #    binCenter.append(hist.GetBinCenter(i+1))
    #    m_eTru.append( hist.GetBinContent(i+1))

    #plt.plot(binCenter, m_eTru, "-")
    plt.plot(binCenter, m_eVis, "--")
    plt.plot(databince, databinco, "o", ms=3)

    plt.show()

        

def fineBinHistPrepare():
    fineBinHist = ROOT.TH1D("fineBin", "", 5000, 10, 60)
    fineBin = np.arange(10, 70, 0.01)
    fineBinContent = []
    for i in range(len(fineBin)):
        fineBinHist.SetBinContent(i-1, theory(fineBin[i], 2000))
        fineBinContent.append( theory(fineBin[i], 2000) )

    #plt.plot(fineBin, fineBinContent, "-" )
    #plt.show()

    return fineBinHist


def main():

    #Ee = np.arange(1, 55, 0.01)
    #pdf = theory1(Ee)
    #pred_npe = []
    #for i in Ee:
    #    pred_npe.append(applyCerNPE(15.7, i))
    #plt.plot(Ee, pdf, "-")
    #plt.show()

    #ke = up.open("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/water/michel/michel_totpe_water1.root")["michel"]["totpe"].array()

    #plt.hist(ke, bins=60, range=(0, 1200))
    #plt.plot(pred_npe, pdf, "-")
    #plt.show()
    #


    #ff = ROOT.TFile("/junofs/users/miaoyu/energy_model/production/J19v1r0-Pre4/water/michel/michel.root", "read")
    #hmc = ff.Get("kine")
    #cent, cont = [], []
    #for i in range(hmc.GetNbinsX()):
    #    cent.append(hmc.GetBinCenter(i))
    #    cont.append(hmc.GetBinContent(i))

    #plt.plot(Ee, pdf, "b-")
    #plt.plot(cent, cont, "r--")
    #plt.show()


    #hist = fineBinHistPrepare()
    #applyCerNPE(hist)

    #applyCerNPE()

    genTheoHist()


if __name__ == "__main__":
    main()
