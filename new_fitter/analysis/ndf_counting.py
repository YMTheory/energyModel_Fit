import ROOT

def countBins(filename, histname):
    ff = ROOT.TFile(filename, "read")
    hist = ff.Get(histname)

    count = 0
    for i in range(hist.GetNbinsX()):
        if hist.GetBinContent(i) != 0:
            count += 1

    return count


def main():

    sources = ["Cs137", "Mn54", "Ge68", "K40", "Co60", "nH", "AmBe", "nC12", "AmC"]

    count = 0 
    for i in sources:
        count += countBins("../output/NewgamB12_kSimQ_kNewAnaCer_kNew/"+i+"hist.root", i+"_data")

    print(count)
    #count += countBins("../output/NewgamB12_kSimQ_kNewAnaCer_kNew/michel.root", "hEvisData")
    #print(count)

    count += countBins("../output/NewgamB12_kSimQ_kNewAnaCer_kNew/spectrum.root", "hData")
    print(count)

    #count += countBins("../spectrum.root", "hData")

    print("Total non-empty bin number = %d" %count)
    


if __name__ == "__main__":
    main()
