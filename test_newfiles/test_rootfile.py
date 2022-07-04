from ROOT import TH2F, TFile, TTree, TCanvas, TH1F, TString
import ROOT
from array import array

ROOT.EnableThreadSafety()
ROOT.EnableImplicitMT()

nHists  = 5 
step    = 20
maxdimx = 1500
maxdimy = 1500
fHist     = []
fCanvas   = []
myFile = TFile.Open("/home/daq/hdd8TB/RSD2/stats_N_script/Run200_W3_croci1.3mm.root", "OPEN")
tree = myFile.Get("Analysis")

nEntries = tree.GetEntries()
for j in range(0, nHists):
	fHist.append(TH2F(str(j),"title",int(maxdimx/step),0,maxdimx,int(maxdimy/step),0,maxdimy))
print("created all histograms")

count = 1
for i in range(0, nEntries):
	if( i == int(nEntries/10*count) ):
		print("done with %d pc of the data"% (count*10))
		count +=1
	tree.GetEntry(i)
	for j in range(0, nHists): fHist[j].Fill(tree.XPos, tree.YPos, tree.ampl[j])
	
for j in range(0, nHists):
	fCanvas.append(TCanvas(str("c%d"%j),str(j),600,600))
	fHist[j].Draw("colz")
	fCanvas[j].Update()
	fCanvas[j].SaveAs(str("c%d.pdf"%j))
