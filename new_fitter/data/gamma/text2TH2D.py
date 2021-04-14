#!/usr/bin/env python
# coding=utf-8

import ROOT
import sys

hh1 = ROOT.TH2D(sys.argv[2], "", 5000, 0, 5000, 100, 0, 100)
hh2 = ROOT.TH2D(sys.argv[3], "", 5000, 0, 5000, 100, 0, 100)

evtid = 0
parid = 0
with open(sys.argv[1]) as f:
    for lines in f.readlines():
        evtid += 1
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
                hh2.SetBinContent(evtid+1, counta, float(j))
            elif "b" in i:
                countb+=1
                tmp = list(i)
                tmp.pop()
                j = ''.join(tmp)
                hh1.SetBinContent(evtid+1, countb, float(j))

ff = ROOT.TFile(sys.argv[4], "recreate")
hh1.Write()
hh2.Write()
ff.Close()

    