import ROOT
import array,numpy,sys
from matplotlib import pyplot as plt

noSys = False

fi = ROOT.TFile.Open("workspace.root")
ws = fi.Get("w")

nbins = 90

rMinimum = -5
rMaximum = 5
np = 100

if noSys: 
  allN = ws.genobj("nuisances")
  iterN = allN.createIterator()
  nNuis = allN.getSize()
  for i in range(nNuis):
    obj = iterN.Next()
    ws.var(obj.GetName()).setConstant(True)

mu   = ws.var("mu"); mu.setConstant()
data = ws.data("obsdata")
pdf  = ws.pdf("combined_pdf")
con  = ws.pdf("nuisance_pdf")

mu.setRange(rMinimum,rMaximum)

data.Print()
con.Print()
pdf.Print()

nll = pdf.createNLL(data,ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(con)))
mu.setConstant(False)
minimF = ROOT.RooMinimizer(nll)
minimF.minimize("Minuit","minimize")
nll0 = nll.getVal()
rmin = mu.getVal()

mu.setConstant(True)
minim = ROOT.RooMinimizer(nll)
# do a scan

myR = [float(rMinimum)+i*float(rMaximum-rMinimum)/np for i in range(np)]
C = []
B = []

for r in myR:
  mu.setVal(r)
  minim.minimize("Minuit","minimize")
  C.append(nll.getVal()-nll0)
  thisB = []
  for b in range(nbins): thisB.append(ws.function("expected_background_bin%d"%(b+1)).getVal())
  B.append(thisB)


# also make the usual tree 
fout = ROOT.TFile("full-LH.root","RECREATE")
tree = ROOT.TTree("limit","limit")

dnll = array.array('d',[0])
r    = array.array('d',[0])
BS = []
for i in range(nbins): BS.append(array.array('d',[0]))

tree.Branch("r",r,"r/D")
tree.Branch("deltaNLL",dnll,"deltaNLL/D")
for i in range(nbins) : tree.Branch("bkg_bin_%d"%(i+1),BS[i],"bkg_bin_%d/D"%(i+1))

#nll0 = min(C)
for i in range(len(myR)) :
  dnll[0]= C[i]
  r[0]=myR[i]
  print r[0],dnll[0]
  for b in range(nbins): BS[b][0]=B[i][b]
  tree.Fill()

fout.cd()
tree.Write()
  

plt.plot(myR,C)
plt.ylabel("-Log(L)")
plt.xlabel("$\mu$")
plt.show()
