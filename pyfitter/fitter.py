from iminuit import cost, Minuit
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

from gammaResponse import gamma
from parameters import parameters


class fitter(object):

    Cs137 = None
    Mn54  = None
    Ge68  = None
    K40   = None
    nH    = None
    Co60  = None
    AmBe  = None
    nC12  = None
    AmC   = None

    bestFit = []
    bestFitError = []

    def __init__(self):
        self.Cs137 = gamma("Cs137", 0.662)
        self.Mn54  = gamma("Mn54", 0.835)
        self.Ge68  = gamma("Ge68", 1.022)
        self.K40   = gamma("K40", 1.461)
        self.nH    = gamma("nH", 2.223)
        self.Co60  = gamma("Co60", 2.506)
        self.AmBe  = gamma("AmBe", 4.430)
        self.nC12  = gamma("nC12", 4.940)
        self.AmC   = gamma("AmC", 6.130)


    def Cs137_pdf(self, x, Ysct, kB, kC, a, b, n):
        print(Ysct, kB, kC, a, b, n)
        # Parameter Setting ...
        parameters.set_Ysct(Ysct)
        parameters.set_kB(kB)
        parameters.set_kC(kC)
        parameters.set_res_a(a)
        parameters.set_res_b(b)
        parameters.set_res_n(n)
    
        self.Cs137.prediction()
        mu, sigma = self.Cs137.get_pred_npe_mean(), self.Cs137.get_pred_npe_sigma()
        return norm.pdf(x, mu, sigma)


    def Mn54_pdf(self, x, Ysct, kB, kC, a, b, n):
        # Parameter Setting ...
        parameters.set_Ysct(Ysct)
        parameters.set_kB(kB)
        parameters.set_kC(kC)
        parameters.set_res_a(a)
        parameters.set_res_b(b)
        parameters.set_res_n(n)
    
        self.Mn54.prediction()
        mu, sigma = self.Mn54.get_pred_npe_mean(), self.Mn54.get_pred_npe_sigma()
        return norm.pdf(x, mu, sigma)


    def Ge68_pdf(self, x, Ysct, kB, kC, a, b, n):
        # Parameter Setting ...
        parameters.set_Ysct(Ysct)
        parameters.set_kB(kB)
        parameters.set_kC(kC)
        parameters.set_res_a(a)
        parameters.set_res_b(b)
        parameters.set_res_n(n)
    
        self.Ge68.prediction()
        mu, sigma = self.Ge68.get_pred_npe_mean(), self.Ge68.get_pred_npe_sigma()
        return norm.pdf(x, mu, sigma)

    def K40_pdf(self, x, Ysct, kB, kC, a, b, n):
        # Parameter Setting ...
        parameters.set_Ysct(Ysct)
        parameters.set_kB(kB)
        parameters.set_kC(kC)
        parameters.set_res_a(a)
        parameters.set_res_b(b)
        parameters.set_res_n(n)
    
        self.K40.prediction()
        mu, sigma = self.K40.get_pred_npe_mean(), self.K40.get_pred_npe_sigma()
        return norm.pdf(x, mu, sigma)
    
    def nH_pdf(self, x, Ysct, kB, kC, a, b, n):
        # Parameter Setting ...
        parameters.set_Ysct(Ysct)
        parameters.set_kB(kB)
        parameters.set_kC(kC)
        parameters.set_res_a(a)
        parameters.set_res_b(b)
        parameters.set_res_n(n)
    
        self.nH.prediction()
        mu, sigma = self.nH.get_pred_npe_mean(), self.nH.get_pred_npe_sigma()
        return norm.pdf(x, mu, sigma)

    def Co60_pdf(self, x, Ysct, kB, kC, a, b, n):
        # Parameter Setting ...
        parameters.set_Ysct(Ysct)
        parameters.set_kB(kB)
        parameters.set_kC(kC)
        parameters.set_res_a(a)
        parameters.set_res_b(b)
        parameters.set_res_n(n)
    
        self.Co60.prediction()
        mu, sigma = self.Co60.get_pred_npe_mean(), self.Co60.get_pred_npe_sigma()
        return norm.pdf(x, mu, sigma)

    def AmBe_pdf(self, x, Ysct, kB, kC, a, b, n):
        # Parameter Setting ...
        parameters.set_Ysct(Ysct)
        parameters.set_kB(kB)
        parameters.set_kC(kC)
        parameters.set_res_a(a)
        parameters.set_res_b(b)
        parameters.set_res_n(n)
    
        self.AmBe.prediction()
        mu, sigma = self.AmBe.get_pred_npe_mean(), self.AmBe.get_pred_npe_sigma()
        return norm.pdf(x, mu, sigma)

    def nC12_pdf(self, x, Ysct, kB, kC, a, b, n):
        # Parameter Setting ...
        parameters.set_Ysct(Ysct)
        parameters.set_kB(kB)
        parameters.set_kC(kC)
        parameters.set_res_a(a)
        parameters.set_res_b(b)
        parameters.set_res_n(n)
    
        self.nC12.prediction()
        mu, sigma = self.nC12.get_pred_npe_mean(), self.nC12.get_pred_npe_sigma()
        return norm.pdf(x, mu, sigma)

    def AmC_pdf(self, x, Ysct, kB, kC, a, b, n):
        # Parameter Setting ...
        parameters.set_Ysct(Ysct)
        parameters.set_kB(kB)
        parameters.set_kC(kC)
        parameters.set_res_a(a)
        parameters.set_res_b(b)
        parameters.set_res_n(n)
    
        self.AmC.prediction()
        mu, sigma = self.AmC.get_pred_npe_mean(), self.AmC.get_pred_npe_sigma()
        return norm.pdf(x, mu, sigma)


    
    def fit(self):
       
        Cs137_data = self.Cs137.get_npe_data()
        Mn54_data = self.Mn54.get_npe_data()
        Ge68_data = self.Ge68.get_npe_data()
        K40_data = self.K40.get_npe_data()
        nH_data = self.nH.get_npe_data()
        Co60_data = self.Co60.get_npe_data()
        AmBe_data = self.AmBe.get_npe_data()
        nC12_data = self.nC12.get_npe_data()
        AmC_data = self.AmC.get_npe_data()

        lf = cost.UnbinnedNLL(Cs137_data, self.Cs137_pdf) + cost.UnbinnedNLL(AmC_data, self.AmC_pdf)
        #lf = cost.UnbinnedNLL(Cs137_data, self.Cs137_pdf) + cost.UnbinnedNLL(Mn54_data, self.Mn54_pdf) + cost.UnbinnedNLL(Ge68_data, self.Ge68_pdf)   + cost.UnbinnedNLL(K40_data, self.K40_pdf) + cost.UnbinnedNLL(nH_data, self.nH_pdf) + cost.UnbinnedNLL(Co60_data, self.Co60_pdf) + cost.UnbinnedNLL(AmBe_data, self.AmBe_pdf) + cost.UnbinnedNLL(nC12_data, self.nC12_pdf) + cost.UnbinnedNLL(AmC_data, self.AmC_pdf)


        m = Minuit(lf, Ysct=1400, kB=6.5e-3, kC=1, a=0.98, b=0.044, n=1.3)
        m.limits["Ysct"] = (1300, 1500)
        m.limits["kB"] = (5.1e-3, 7.5e-3)
        m.limits["kC"] = (0.5, 1.5)
        m.limits["a"] = (0.5, 1.5)
        m.limits["b"] = (0.0, 0.1)
        m.limits["n"] = (1.0, 2.0)
        m.migrad()
        m.hesse()
        
        print("===== Fitting Results =====")
        print(m.values)
        print(m.errors)
        print(m.covariance)

        for i, j in zip(m.values, m.errors):
            self.bestFit.append(i)
            self.bestFitError.append(j)

    def plot(self):
        Cs137_data = self.Cs137.get_npe_data()
        Mn54_data = self.Mn54.get_npe_data()
        Ge68_data = self.Ge68.get_npe_data()
        K40_data = self.K40.get_npe_data()
        nH_data = self.nH.get_npe_data()
        Co60_data = self.Co60.get_npe_data()
        AmBe_data = self.AmBe.get_npe_data()
        nC12_data = self.nC12.get_npe_data()
        AmC_data = self.AmC.get_npe_data()

        fig = plt.figure(figsize=(12, 10))
        spec = gridspec.GridSpec(ncols=3, nrows=3)

        ax0 = fig.add_subplot(spec[0])
        ax1 = fig.add_subplot(spec[1])
        ax2 = fig.add_subplot(spec[2])
        ax3 = fig.add_subplot(spec[3])
        ax4 = fig.add_subplot(spec[4])
        ax5 = fig.add_subplot(spec[5])
        ax6 = fig.add_subplot(spec[6])
        ax7 = fig.add_subplot(spec[7])
        ax8 = fig.add_subplot(spec[8])

        ax0.hist(Cs137_data, bins=100, density=True, histtype="step")
        x = np.arange(np.min(Cs137_data), np.max(Cs137_data), 1)
        ax0.plot(x, self.Cs137_pdf(x, *self.bestFit), "-")
        ax0.set_xlabel("NPE", fontsize=13)
        ax0.set_ylabel("count", fontsize=13)

        ax1.hist(Mn54_data, bins=100, density=True, histtype="step")
        x = np.arange(np.min(Mn54_data), np.max(Mn54_data), 1)
        ax1.plot(x, self.Mn54_pdf(x, *self.bestFit), "-")
        ax0.set_xlabel("NPE", fontsize=13)
        ax0.set_ylabel("count", fontsize=13)

        ax2.hist(Ge68_data, bins=100, density=True, histtype="step")
        x = np.arange(np.min(Ge68_data), np.max(Ge68_data), 1)
        ax2.plot(x, self.Ge68_pdf(x, *self.bestFit), "-")
        ax2.set_xlabel("NPE", fontsize=13)
        ax2.set_ylabel("count", fontsize=13)

        ax3.hist(K40_data, bins=100, density=True, histtype="step")
        x = np.arange(np.min(K40_data), np.max(K40_data), 1)
        ax3.plot(x, self.K40_pdf(x, *self.bestFit), "-")
        ax3.set_xlabel("NPE", fontsize=13)
        ax3.set_ylabel("count", fontsize=13)

        ax4.hist(nH_data, bins=100, density=True, histtype="step")
        x = np.arange(np.min(nH_data), np.max(nH_data), 1)
        ax4.plot(x, self.nH_pdf(x, *self.bestFit), "-")
        ax4.set_xlabel("NPE", fontsize=13)
        ax4.set_ylabel("count", fontsize=13)

        ax5.hist(Co60_data, bins=100, density=True, histtype="step")
        x = np.arange(np.min(Co60_data), np.max(Co60_data), 1)
        ax5.plot(x, self.Co60_pdf(x, *self.bestFit), "-")
        ax5.set_xlabel("NPE", fontsize=13)
        ax5.set_ylabel("count", fontsize=13)

        ax6.hist(AmBe_data, bins=100, density=True, histtype="step")
        x = np.arange(np.min(AmBe_data), np.max(AmBe_data), 1)
        ax6.plot(x, self.AmBe_pdf(x, *self.bestFit), "-")
        ax6.set_xlabel("NPE", fontsize=13)
        ax6.set_ylabel("count", fontsize=13)

        ax7.hist(nC12_data, bins=100, density=True, histtype="step")
        x = np.arange(np.min(nC12_data), np.max(nC12_data), 1)
        ax7.plot(x, self.nC12_pdf(x, *self.bestFit), "-")
        ax7.set_xlabel("NPE", fontsize=13)
        ax7.set_ylabel("count", fontsize=13)

        ax8.hist(AmC_data, bins=100, density=True, histtype="step")
        x = np.arange(np.min(AmC_data), np.max(AmC_data), 1)
        ax8.plot(x, self.AmC_pdf(x, *self.bestFit), "-")
        ax8.set_xlabel("NPE", fontsize=13)
        ax8.set_ylabel("count", fontsize=13)

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=0.02, hspace=0.02)
        plt.tight_layout()

        plt.savefig("gamma_fit.pdf")
    


