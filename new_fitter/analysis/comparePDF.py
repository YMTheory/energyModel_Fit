import prmBetaLoader as bloader
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import ROOT


def add_subplot_axes(ax, rect, axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    #subax = fig.add_axes([x,y,width,height],facecolor=facecolor)  # matplotlib 2.0+
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax




def main():

    ##### Livermore
    filename1 = "../data/gamma/log-livermore.txt"
    print(filename1)
    prmBeta1, prmAntiBeta1 = bloader.loadPrmBeta(filename1)

    allBeta1 = []
    BetaNum1 = np.zeros(100000)
    for i in range(100000):
        for j in prmBeta1[i]:
            allBeta1.append(j)
            BetaNum1[i] += 1
        for k in prmAntiBeta1[i]:
            allBeta1.append(k)
            BetaNum1[i] += 1
    
    

    ##### Penelope
    filename2 = "../data/gamma/log-penelope.txt"
    print(filename2)
    prmBeta2, prmAntiBeta2 = bloader.loadPrmBeta(filename2)

    allBeta2 = []
    BetaNum2 = np.zeros(100000)
    for i in range(100000):
        for j in prmBeta2[i]:
            allBeta2.append(j)
            BetaNum2[i] += 1
        for k in prmAntiBeta2[i]:
            allBeta2.append(k)
            BetaNum2[i] += 1
    

    """
    ### chi2 TEST
    h1 = ROOT.TH1D("h1", "h1", 100, 0, 8)
    h2 = ROOT.TH1D("h2", "h2", 100, 0, 8)
    for i, j in zip(allBeta1, allBeta2):
        h1.Fill(i)
        h2.Fill(j)
    chi2 = h1.KolmogorovTest(h2)
    print("EnergyPdf chi2test : ", chi2)
    
    h3 = ROOT.TH1I("h3", "h3", 50, 0, 50)
    h4 = ROOT.TH1I("h4", "h4", 50, 0, 50)
    for i, j in zip(BetaNum1, BetaNum2):
        h3.Fill(i)
        h4.Fill(j)
    chi2 = h3.KolmogorovTest(h4)
    print("NumberPdf chi2test : ", chi2)

    """



    #fig = plt.figure(figsize=(6, 9))
    #spec = gridspec.GridSpec(ncols=1, nrows=4)


    #fig = plt.figure(constrained_layout=True, figsize=(6, 10))
    #heights = [1, 2, 1, 2]
    #spec = fig.add_gridspec(ncols=1, nrows=4, height_ratios=heights)
    #ax0 = fig.add_subplot(spec[1])
    #ax1 = fig.add_subplot(spec[3])
    #ax2 = fig.add_subplot(spec[0])
    #ax3 = fig.add_subplot(spec[2])
    fig = plt.figure(constrained_layout=True, figsize=(10, 6))
    heights = [1, 2]
    spec = fig.add_gridspec(ncols=2, nrows=2, height_ratios=heights)
    ax0 = fig.add_subplot(spec[3])
    ax2 = fig.add_subplot(spec[1])
    ax1 = fig.add_subplot(spec[2])
    ax3 = fig.add_subplot(spec[0])



    cont1, edge1 = np.histogram(allBeta1, bins=100, range=(0, 8))
    cont2, edge2 = np.histogram(allBeta2, bins=100, range=(0, 8))
    bin1 = np.linspace(0, 8, 100)


    ax0.hist(allBeta1, bins=100, histtype="step", color="blue", label="Livermore")
    ax0.hist(allBeta2, bins=100, histtype="step", color="red", label="Penelope")

    #ax0.bar(bin1, cont1/np.sum(cont1), width=8/100., color="white", edgecolor="blue", label="Livermore")
    #ax0.bar(bin1, cont2/np.sum(cont2), width=8/100., color="white", edgecolor="red", label="Penelope")

    ax0.set_xlabel("(b) Primary beta energy [MeV]", fontsize=15)
    ax0.set_ylabel("counts per bin", fontsize=15)
    ax0.legend(prop={"size":14}, ncol=2)
    ax0.semilogy()
    ax0.set_ylim(0, 5000000)
    ax0.tick_params(axis='both', which='major', labelsize=14)
    ax0.set_xlim(-0.1, 8)


    #rect1 = [0.35, 0.52, 0.6, 0.4]
    #ax2 = add_subplot_axes(ax0, rect1)
    ax2.errorbar(bin1, (cont2-cont1)/cont2, yerr=np.sqrt(cont1/cont2**2+cont2*cont1**2/cont2**4), fmt="o", color="black", ms=3, mfc="w")
    #ax2.set_xlabel("Primary beta energy [MeV]", fontsize=13)
    ax2.set_ylabel("Bias", fontsize=13)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.grid(True)
    ax2.set_ylim(-0.3, 0.3)
    ax2.set_xlim(-0.1, 8)


    ax1.hist(BetaNum1, bins=50, range=(0, 50), histtype="step", color="blue", label="Livermore")
    ax1.hist(BetaNum2, bins=50, range=(0, 50), histtype="step", color="red", label="Penelope")

    ax1.set_xlabel("(a) Primary beta multiplicity", fontsize=15)
    ax1.set_ylabel("counts per bin", fontsize=15)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.legend(loc="lower center", prop={"size":14}, ncol=1)
    ax1.semilogy()
    ax1.set_ylim(0, 20000)
    ax1.set_xlim(0, 50)



    #rect2 = [0.21, 0.24, 0.5, 0.33]
    #ax3 = add_subplot_axes(ax1, rect2)
    cont1, edge1 = np.histogram(BetaNum1, bins=50, range=(0, 50))
    cont2, edge2 = np.histogram(BetaNum2, bins=50, range=(0, 50))
    bin1 = np.linspace(0, 50, 50)
    ax3.errorbar(bin1, (cont2-cont1)/cont2, yerr=np.sqrt(cont1/cont2**2+cont2*cont1**2/cont2**4), fmt="o", color="black", ms=3, mfc="w")
    #ax3.set_xlabel("Primary beta multiplicity", fontsize=13)
    ax3.set_ylabel("Bias", fontsize=13)
    ax3.tick_params(axis='both', which='major', labelsize=12)
    ax3.grid(True)
    ax3.set_ylim(-0.3, 0.3)
    ax3.set_xlim(0, 50)


    plt.subplots_adjust(left=None, bottom=None, right=None, top=None , wspace=None, hspace=None)

    plt.tight_layout()
    plt.savefig("compareModel.pdf")
    plt.show()



if __name__ == "__main__" :
    main()
