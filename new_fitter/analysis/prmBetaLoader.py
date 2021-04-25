
def loadPrmBeta(filename):
    secBetaArr, secAntiBetaArr = [], []
    lineId = 0
    with open(filename) as f:
        for lines in f.readlines():
            #if lineId % 50 == 0:
            #    print("Process %d percent ..." %int(lineId/50))
            oneEvtBeta, oneEvtAntiBeta = [], []
            line = lines.strip("\n")
            data = line.split(" ")
            counta = 0
            countb = 0
            for i in data:
                if "a" in i:
                    counta+=1
                    tmp = list(i)
                    tmp.pop()
                    j = ''.join(tmp)
                    oneEvtAntiBeta.append(float(j))
                #    hh2.SetBinContent(evtid+1, counta, float(j))
                if "b" in i:
                    countb+=1
                    tmp = list(i)
                    tmp.pop()
                    j = ''.join(tmp)
                    oneEvtBeta.append(float(j))
            #secBetaArr.append(tmpE)
            lineId += 1
            secAntiBetaArr.append(oneEvtAntiBeta)
            secBetaArr.append(oneEvtBeta)
    return secBetaArr, secAntiBetaArr


