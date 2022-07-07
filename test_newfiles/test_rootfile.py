import sys
import ROOT
from ROOT import TH2F, TFile, TTree, TCanvas, TH1F, TString
from array import array
from progress.bar import ChargingBar

ROOT.EnableThreadSafety()
ROOT.EnableImplicitMT(8)

nHists  = 16 
step    = 10
maxdimx = 1100
maxdimy = 700
fHist     = []
fCanvas   = []

if len(sys.argv) == 1:
	filename = "../Run22_150V.root" #"../RunXX.root" # "/home/daq/hdd8TB/RSD2/stats_N_script/Run22_240.root" # " 
else:
	filename = sys.argv[1]
myFile = TFile.Open(filename, "OPEN")
tree = myFile.Get("Analysis")

nEntries = tree.GetEntries()
for j in range(nHists):
	fHist.append(TH2F(str(j),"c"+str(j),int(maxdimx/step),0,maxdimx,int(maxdimy/step),0,maxdimy))
print("created all histograms")

bar = ChargingBar('Processing', max=nEntries, suffix = '%(percent)d%% [%(elapsed_td)s]')
for i in range(nEntries):
	tree.GetEntry(i)
	for j in range(nHists): fHist[j].Fill(tree.XPos, tree.YPos, tree.ampl[j])
	bar.next()
bar.finish()

for j in range(nHists):
	fCanvas.append(TCanvas(str("c%d"%j),str(j),600,600))
	fHist[j].Draw("colz")
	fCanvas[j].Modified()
	fCanvas[j].Update()
	fCanvas[j].SaveAs(str("c%d.pdf"%j))
