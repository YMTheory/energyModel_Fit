import numpy as np

class truthLoader(object):

    def __init__(self, dirs):
        self.parTruth = {}
        self.det = dirs

    def setDet(self, det):
        self.det = det

    def setDet(self, det):
        self.det = det

    def getNcer(self, mom):
        keyname = str(mom) 
        return self.parTruth[keyname]["Ncer"]

    def getNPE(self, mom):
        keyname = str(mom) 
        return self.parTruth[keyname]["NPE"]

    def getNPEerr(self, mom):
        keyname = str(mom) 
        return self.parTruth[keyname]["NPEerr"]

    def getsigma(self, mom):
        keyname = str(mom) 
        return self.parTruth[keyname]["sigma"]

    def getsigmaerr(self, mom):
        keyname = str(mom) 
        return self.parTruth[keyname]["sigmaerr"]


    def fitElecNPESimTruth(self, KE):
        import uproot as up
        import ROOT
        for mom in KE:
            filename = "/junofs/users/miaoyu/simulation/LS_Sim/jobs/"+self.det+"/e-/newsim" + str(int(mom*1000)) + ".root"
            f = up.open(filename)
            cerpe = f["photon"]["cerPE"].array()
            totpe = f["photon"]["totPE"].array()

            h = ROOT.TH1D("h", "", 100, np.min(totpe), np.max(totpe))
            for n in totpe:
                h.Fill(n)
            h.Fit("gaus", "Q0")
            mu = h.GetFunction("gaus").GetParameter(1)
            muerr = h.GetFunction("gaus").GetParError(1)
            sig = h.GetFunction("gaus").GetParameter(2)
            sigerr = h.GetFunction("gaus").GetParError(2)

            tmp_dict = {"Ncer" : 0}
            tmp_dict["Ncer"] = np.mean(cerpe)
            tmp_dict["NPE"] = mu
            tmp_dict["NPEerr"] = muerr
            tmp_dict["sigma"] = sig
            tmp_dict["sigmaerr"] = sigerr

            keyname = str(mom) + "MeV"
            self.parTruth[keyname] = tmp_dict

            del h



    def fitPosiNPESimTruth(self, KE):
        import uproot as up
        import ROOT
        for mom in KE:
            filename = "/junofs/users/miaoyu/simulation/LS_Sim/jobs/"+self.det+"/e+/newsim" + str(int(mom*1000)) + ".root"
            f = up.open(filename)
            cerpe = f["photon"]["cerPE"].array()
            totpe = f["photon"]["totPE"].array()

            h = ROOT.TH1D("h", "", 100, np.min(totpe), np.max(totpe))
            for n in totpe:
                h.Fill(n)
            h.Fit("gaus", "Q0")
            mu = h.GetFunction("gaus").GetParameter(1)
            muerr = h.GetFunction("gaus").GetParError(1)
            sig = h.GetFunction("gaus").GetParameter(2)
            sigerr = h.GetFunction("gaus").GetParError(2)

            tmp_dict = {"Ncer" : 0}
            tmp_dict["Ncer"] = np.mean(cerpe)
            tmp_dict["NPE"] = mu
            tmp_dict["NPEerr"] = muerr
            tmp_dict["sigma"] = sig
            tmp_dict["sigmaerr"] = sigerr

            keyname = str(mom) + "MeV"
            self.parTruth[keyname] = tmp_dict

            del h



    def fitGammaNPESimTruth(self):
        import uproot as up
        import ROOT
        for name in ["Cs137", "Mn54", "K40", "nH", "AmBe", "AmC"]:
            filename = "/junofs/users/miaoyu/simulation/LS_Sim/jobs/"+self.det+"/gamma/" + name  + "_new.root"
            f = up.open(filename)
            totpe = f["photon"]["totPE"].array()

            h = ROOT.TH1D("h", "", 100, np.min(totpe), np.max(totpe))
            for n in totpe:
                h.Fill(n)
            h.Fit("gaus", "Q0")
            mu = h.GetFunction("gaus").GetParameter(1)
            muerr = h.GetFunction("gaus").GetParError(1)
            sig = h.GetFunction("gaus").GetParameter(2)
            sigerr = h.GetFunction("gaus").GetParError(2)

            tmp_dict = {"Ncer" : 0}
            tmp_dict["Ncer"] = 0
            tmp_dict["NPE"] = mu
            tmp_dict["NPEerr"] = muerr
            tmp_dict["sigma"] = sig
            tmp_dict["sigmaerr"] = sigerr

            keyname = name
            self.parTruth[keyname] = tmp_dict

            del h

