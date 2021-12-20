import numpy as np
import matplotlib.pyplot as plt
import elecLoader as eloader
from matplotlib import gridspec
import ROOT

def func(p0, p1, p2, x):
    s2 = p0 + p1*x + p2 * x**2
    return s2 if s2>0  else 0
    #if s2 < 0:
    #    return 0
    #else:
    ##    #return np.sqrt(s2)
    #    return s2




def loadCerResFile(filename):
    E, cerpe, sigma = [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            E.append(float(data[0]))
            cerpe.append(float(data[1]))
            sigma.append(float(data[2]))

    E = np.array(E)
    cerpe = np.array(cerpe)
    sigma = np.array(sigma)

    return E, cerpe, sigma
            
def loadCovFile(filename):
    E, cerpe, sigma, sigma_err = [], [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            E.append(float(data[0]))
            cerpe.append(float(data[1]))
            sigma.append(float(data[2]))
            sigma_err.append(float(data[3]))

    E = np.array(E)
    cerpe = np.array(cerpe)
    sigma = np.array(sigma)
    sigma_err = np.array(sigma_err)

    return E, cerpe, sigma, sigma_err


def sample_uncorelation(filename, num, sampleSize):

    num_sample = sampleSize

    cov = loadCov(filename, num)

    pars = [ np.zeros((1, sampleSize)) for i in range(num) ]
    for i in range(num):
        pars[i] = np.random.normal(loc=0, scale=np.sqrt(scale=cov[i ,i]), size=sampleSize )

    return pars
        


from scipy.linalg import eigh, cholesky
from scipy.stats import norm

def loadBestFit(filename, num):
    bestFit, err = [], []
    row = 0
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if row == num:
                for kk in range(num):
                    bestFit.append(float(data[kk]))
            if row == num+1:
                for kk in range(num):
                    err.append(float(data[kk]))
            row += 1

    return bestFit, err


def loadCov(filename, num):
    row, col = 0, 0
    cov_mat = np.ones((num, num))
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if row == num:
                break
            col = 0
            for i in data:
               cov_mat[row, col] = float(i) 
               cov_mat[col, row] = float(i)
               col+=1
            row += 1

    return cov_mat


def sample_corelation(filename, num, sampleSize):
    method = 'eigenvectors'
    
    num_sample = sampleSize

    cov = loadCov(filename, num)
    #print(cov)
    
    x = norm.rvs(size=(num, num_sample))

    if method == 'cholesky':
        # Compute the Cholesky decomposition.
        c = cholesky(cov, lower=True)
    else:
        # Compute the eigenvalues and eigenvectors.
        evals, evecs = eigh(cov)
        # Construct c, so c*c^T = r.
        c = np.dot(evecs, np.diag(np.sqrt(evals)))

    y = np.dot(c, x)

    return y    


def resFunc(x, a, b, c):
    s2 = a + b*x + c*x**2
    if s2<0:
        return 0
    else:
        return np.sqrt(s2)





def main():

    es = 3134.078/2.223

    draw = []
    for i in np.arange(100, 9200, 100):
        draw.append(i)
    draw = np.array(draw)
    
    
    #a1, b1, c1 = 9.89870e-01, 7.78072e-03, 0
    #a1, b1, c1 = 9.9042e-1, 7.91898e-03, 0
    #a1, b1, c1 = 0.988, 7.89e-3, 0
    #a1err, b1err = 6.46e-3, 4.33e-4
    #filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12_kSimQ_kNewAnaC1_kNPE/gamB12_kSimulation_kNewAnaC1_kNPE_rescov.txt"
    #filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12Mic_kSimQ_kNewAnaC1_kNPE/gamB12Mic_kSimulation_kNewAnaCer_kNPE_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12Mic_kSimQ_kNewAnaC1_kNPE/gamB12Mic_kSimulation_kNewAnaCer_kNPE_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_gamB12Mic_kSimQ_kNewAnaC1_kNPE_2/gamB12Mic_kSimulation_kNewAnaCer_kNPE_rescov2.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/articles_Mic_kSimQ_kNewAnaC1_kNPE/Mic_kSimulation_kNewAnaCer_kNPE_rescov3.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/gamB12NewMicwhole/gamB12NewMicwhole_kSimQ_kCerNewAna_kNPE_rescov.txt"
    #filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/gamB12NewMicedge/gamB12NewMicedge_kSimQ_kCerNewAna_kNPE_rescov.txt"
    #filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/gamB12NewMic1whole/gamB12NewMic1whole_kSimQ_kCerNewAna_kNPE_rescov.txt"
    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/Newgam/Newgam_kSimQ_kNewCerAna_kNPE_rescov.txt"
    filename3 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12NewMic1whole/NewgamB12NewMicwhole1_rescov.txt"
    filename2 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamB12/NewgamB12_kSimQ_kNewCerAna_kNPE_rescov.txt"
    filename1 = "/junofs/users/miaoyu/energy_model/energyModel_Fit/new_fitter/output/NewgamNewB12_kSimQ_kNewAnaCer_kNew/NewgamNewB12_kSimQ_kNewAnaCer_kNew_rescov.txt"


    par, parerr = loadBestFit(filename1, 2)
    a1, b1 = par[0], par[1]
    c1 = 0
    a1err, b1err = parerr[0], parerr[1]
    p1 = a1*a1
    p2 = b1*b1
    p1err, p2err = 2*a1*a1err, 2*b1*b1err

    par, parerr = loadBestFit(filename2, 2)
    a2, b2 = par[0], par[1]
    c2 = 0
    a2err, b2err = parerr[0], parerr[1]
    p3 = a2*a2
    p4 = b2*b2
    p3err, p4err = 2*a2*a2err, 2*b2*b2err

    par, parerr = loadBestFit(filename3, 2)
    a3, b3 = par[0], par[1]
    c3 = 0
    a3err, b3err = parerr[0], parerr[1]
    p5 = a3*a3
    p6 = b3*b3
    p5err, p6err = 2*a3*a3err, 2*b3*b3err


    ymin1, ymax1 = [], []
    ymin2, ymax2 = [], []
    ymin3, ymax3 = [], []
    sampleSize = 5000
    sigma1 = sample_corelation(filename1, 2, sampleSize)
    sigma2 = sample_corelation(filename2, 2, sampleSize)
    sigma3 = sample_corelation(filename3, 2, sampleSize)
    dnum = len(draw)
    for i in range(dnum):
        ymin1.append(1000000)
        ymax1.append(-100)
        ymin2.append(1000000)
        ymax2.append(-100)
        ymin3.append(1000000)
        ymax3.append(-100)
    
    for i in range(sampleSize):
    
        m_a1 = a1 + sigma1[0, i]
        m_b1 = b1 + sigma1[1, i]
        m_c1 = 0
        m_a2 = a2 + sigma2[0, i]
        m_b2 = b2 + sigma2[1, i]
        m_c2 = 0
        m_a3 = a3 + sigma3[0, i]
        m_b3 = b3 + sigma3[1, i]
        m_c3 = 0

        #m_a3 = ROOT.gRandom.Gaus(a1, a1err)
        #m_b3 = ROOT.gRandom.Gaus(b1, b1err)
        #m_c3 = 0

        totsigma1 = []
        totsigma2 = []
        totsigma3 = []
        for i in draw:
            totsigma1.append( np.sqrt(func(0, m_a1*m_a1, m_b1*m_b1, i) ) )
            totsigma2.append( np.sqrt(func(0, m_a2*m_a2, m_b2*m_b2, i) ) )
            totsigma3.append( np.sqrt(func(0, m_a3*m_a3, m_b3*m_b3, i) ) )
        for n in range(dnum):
            if totsigma1[n] < ymin1[n]:
                ymin1[n] = totsigma1[n]
            if totsigma1[n] > ymax1[n]:
                ymax1[n] = totsigma1[n]
            if totsigma2[n] < ymin2[n]:
                ymin2[n] = totsigma2[n]
            if totsigma2[n] > ymax2[n]:
                ymax2[n] = totsigma2[n]
            if totsigma3[n] < ymin3[n]:
                ymin3[n] = totsigma3[n]
            if totsigma3[n] > ymax3[n]:
                ymax3[n] = totsigma3[n]
    
    
    ymin1 = np.array(ymin1)
    ymax1 = np.array(ymax1)
    ymin2 = np.array(ymin2)
    ymax2 = np.array(ymax2)
    ymin3 = np.array(ymin3)
    ymax3 = np.array(ymax3)
    
    


    resolE, totpe, resolData, resolerr = eloader.getResolArray()

    resolE1, resolData1, resolerr1 = [], [], []
    with open("../data/electron/elecResol4.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if float(data[0]) > 20:
                continue
            resolE1.append(float(data[1]))
            resolData1.append(float(data[5]))
            resolerr1.append(float(data[6]))


    resolFit  = []
    #resolModel = []
    resolFit1 =  []
    resolFit2 =  []
    resolFit3 = []
    #for i in etrue:
    #    resolFit.append( np.sqrt(func(p0, p1, p2, i)) )
    #    resolModel.append( np.sqrt(func(a0, a1, a2, i)) )
    #    resolFit1.append( np.sqrt(func(-123.213, 3.39232, 2.53173e-3, i)) )
    resolE2 = []
    for i in draw:
        resolE2.append(i/es)
        resolFit1.append( np.sqrt(func(0, p1, p2, i)) )
        resolFit2.append( np.sqrt(func(0, p3, p4, i)) )
        resolFit3.append( np.sqrt(func(0, p5, p6, i)) )



    resolE1 = np.array(resolE1)
    resolData1 = np.array(resolData1)
    resolerr1 = np.array(resolerr1)
    resolFit1 = np.array(resolFit1)
    resolFit2 = np.array(resolFit2)
    resolFit3 = np.array(resolFit3)

    


    #resAlone = []
    #for i in draw:
    #    resAlone.append(np.sqrt(func(0, 1.086**2, 4.321e-03**2, i)))
    
    fig, ax = plt.subplots()
    #fig = plt.figure(figsize=(10, 5))
    #spec = gridspec.GridSpec(ncols=2, nrows=1)
    #                     #height_ratios=[1, 2])

    #ax = fig.add_subplot(spec[0])
    #ax0 = fig.add_subplot(spec[1])

    #ax0.plot(draw/es, resolFit1/resAlone, "-", color="orange", lw=2)
    #ax0.errorbar(resolE1/es, resolData1, yerr=resolerr1, fmt="o", ms=6, lw=2, color='black', mfc="w", label="Simulation truth")
    #ax0.plot(draw/es, resolFit1/draw, "-", color="crimson", zorder=2, lw=2, label="only Michel_whole")
    #ax0.plot(draw/es, resAlone/draw, "-",  color="blue", zorder=2, lw=2, label="Parameterization")
    #ax0.plot(draw/es, resolFit3/draw, "--",  color="green", zorder=2, lw=2, label="gam+Michel_whole+B12")
    #ax0.set_ylabel("Ratio", fontsize=17, color="black")
    #ax0.set_xlim(0, 6)
    #ax0.tick_params(axis='y', which='major', labelsize=15, labelcolor="black")
    #ax0.tick_params(axis='x', which='major', labelsize=0)
    #ax0.grid(True)


    ax.errorbar(resolE1/es, resolData1, yerr=resolerr1, fmt="o", ms=6, lw=2, color='black', mfc='w', label="Simulation truth")
    ax.plot(draw/es, resolFit1/draw, "-", color="crimson", zorder=2, lw=2, label="gamma only")
    ax.plot(draw/es, resolFit2/draw, "--",  color="blue", zorder=2, lw=2, label=r"gamma+B12")
    ax.legend(prop={'size':14})
    ax.fill_between(draw/es, (resolFit1-10*(resolFit1-ymin1))/draw, (resolFit1 + 10*(ymax1-resolFit1))/draw, color="crimson", alpha=0.3)
    ax.fill_between(draw/es, (resolFit2-10*(resolFit2-ymin2))/draw, (resolFit2 + 10*(ymax2-resolFit2))/draw, color="royalblue", alpha=0.3)

    ax.set_xlabel(r"Electron  $E^{vis}$ [MeV]", fontsize=14)
    ax.set_ylabel(r"${\sigma}/{E^{vis}}$", fontsize=15, color="black")
    ax.tick_params(axis='y', which='major', labelsize=13, labelcolor="black")
    ax.tick_params(axis='x', which='major', labelsize=13)
    ax.grid(True)



    #ax0.plot(draw/es, (resolFit3-resolFit1)/resolFit3, "-", lw=2, color="black")
    ##ax0.errorbar(resolE1/es, resolData1, yerr=resolerr1, fmt="o", ms=6, lw=2, color='black', mfc='w', label="Simulation truth")
    ##ax0.plot(draw/es, resolFit1/draw, "-", color="crimson", zorder=2, lw=2, label="gamma only")
    ##ax0.plot(draw/es, resolFit3/draw, "--",  color="black", zorder=2, lw=2, label=r"gamma+B12+Michel")
    ##ax0.legend(prop={'size':14})
    ##ax0.fill_between(draw/es, (resolFit1-10*(resolFit1-ymin1))/draw, (resolFit1 + 10*(ymax1-resolFit1))/draw, color="crimson", alpha=0.3)
    ##ax0.fill_between(draw/es, (resolFit3-10*(resolFit3-ymin3))/draw, (resolFit3 + 10*(ymax3-resolFit3))/draw, color="slategrey", alpha=0.5)

    #ax0.set_xlabel(r"Electron  $E^{vis}$ [MeV]", fontsize=14)
    ##ax0.set_ylabel(r"${\sigma}/{E^{vis}}$", fontsize=15, color="black")
    #ax0.set_ylabel("Diff", fontsize=15, color="black")
    #ax0.tick_params(axis='y', which='major', labelsize=13, labelcolor="black")
    #ax0.tick_params(axis='x', which='major', labelsize=13)
    #ax0.grid(True)

    #textstr = '\n'.join((
    #r'$a = 9.9\times10^{-1} \pm 6.0\times10^{-3}$' ,
    #r'$b = 7.8\times10^{-3} \pm 4.3\times10^{-5}$' 
    #))

    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #ax.text(3, 0.0012, textstr, fontsize=15, bbox=props)

    plt.tight_layout()

    #plt.savefig("sigmaE_elec3.pdf")
    plt.show()




if __name__ == "__main__":
    main()




