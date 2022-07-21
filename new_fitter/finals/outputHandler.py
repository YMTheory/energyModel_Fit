import numpy as np

class outputHandler(object):

    def __init__(self, filename, Npar):
        self.filename = filename
        self.Npar = Npar
        self.sampleSize = 5000

    
    def loadBestFit(self):
        """ load best fit values and errors """
        bestFit, err = [], []
        row = 0
        with open(self.filename) as f:
            for lines in f.readlines():
                line = lines.strip("\n")
                data = line.split(" ")
                if row == self.Npar:
                    for kk in range(self.Npar):
                        bestFit.append(float(data[kk]))
                if row == self.Npar+1:
                    for kk in range(self.Npar):
                        err.append(float(data[kk]))
                row += 1

        return bestFit, err


    def loadCov(self):
        row, col = 0, 0
        cov_mat = np.ones((self.Npar, self.Npar))
        with open(self.filename) as f:
            for lines in f.readlines():
                line = lines.strip("\n")
                data = line.split(" ")
                if row == self.Npar:
                    break
                col = 0
                for i in data:
                   cov_mat[row, col] = float(i) 
                   cov_mat[col, row] = float(i)
                   col+=1
                row += 1
    
        return cov_mat
    
    def sample_corelation(self):
        from scipy.linalg import eigh, cholesky
        from scipy.stats import norm
        method = 'eigenvectors'
        
        num_sample = self.sampleSize
    
        cov = self.loadCov()
        #print(cov)
        
        x = norm.rvs(size=(self.Npar, num_sample))
    
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



    def getSampleSize(self):
        return self.sampleSize



